function hlp_worker(varargin)
% Act as a lightweight worker process for use with hlp_schedule.
% hlp_worker(Options...)
%
% Receives commands (string expressions) from the network, evaluate them, and send off the result to 
% some collector (again as a string). Processing is done in a single thread.
%
% In:
%   Options... : optional name-value pairs, with possible names:
%                 'port': port number on which to listen for requests (default: 23547)
%                         if the port is already in use, the next free one will be chosen,
%                         until port+portrange is exceeded; then, a free one will be chosen
%                         if specified as 0, a free one is chosen directly
%
%                 'portrange': number of ports to try following the default/supplied port (default: 16)
%
%                 'backlog': backlog of queued incoming connections (default: 0)
%
%                 'timeout_accept': timeout for accepting connections, in seconds (default: 3)
%
%                 'timeout_send': timeout for sending results, in seconds (default: 10)
%
%                 'timeout_recv': timeout for receiving data, in seconds (default: 5)
%
% Notes:
%  * use multiple workers to make use of multiple cores
%  * use only ports that are not accessible from the internet
%  * request format: <task_id><collectoraddress_length><collectoraddress><body_length><body>
%    <task_id>: identifier of the task (needs to be forwarded, with the result,
%               to a some data collector upon task completion) (int)
%    <collectoraddress_length>: length, in bytes, of the data collector's address (int)
%    <collectoraddress>: where to send the result for collection (string, formatted as host:port)
%    <body_length>: length, in bytes, of the message body (int)
%    <body>: a MATLAB command that yields, when evaluated, some result in ans (string, as MATLAB expression)
%            if an exception occurs, the exception struct (as from lasterror) replaces ans
%  * response format: <task_id><body_length><body>
%
%                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                       2010-08-26

import java.io.*
import java.net.*

% read options
opts = hlp_varargin2struct(varargin, 'port',23547, 'portrange',16, 'backlog',1, 'timeout_accept',3, ...
    'timeout_send',10, 'timeout_recv',5, 'receive_buffer',64000);

% open a new server socket (first trying the specified portrange, then falling back to 0)
for port = [opts.port:opts.port+opts.portrange 0]
    try
        serv = ServerSocket(port, opts.backlog);
        break;
    catch,end
end
disp(['This is ' hlp_hostname ' (' hlp_hostip '). Listening on port ' num2str(serv.getLocalPort())]);
% set socket properties (e.g., making the function interruptible)
serv.setReceiveBufferSize(opts.receive_buffer);
serv.setSoTimeout(round(1000*opts.timeout_accept));
% make sure that the server socket will be closed when this function is terminated
cleaner = onCleanup(@()serv.close());

tasknum = 1;
disp('waiting for connections...');
while 1
    try
        % wait for an incoming request        
        conn = serv.accept();
        conn.setSoTimeout(round(1000*opts.timeout_recv));
        conn.setTcpNoDelay(1);
        disp('connected.');
                
        try
            % parse request
            in = DataInputStream(conn.getInputStream());
            cr = ChunkReader(in);
            taskid = in.readInt();
            collector = char(cr.readFully(in.readInt())');
            task = char(cr.readFully(in.readInt())');
            disp('received data; replying.');
            out = DataOutputStream(conn.getOutputStream());
            out.writeInt(taskid+length(collector)+length(task));
            out.flush();
            conn.close();
            
            % evaluate task & serialize result
            disp(['running task ' num2str(taskid) ' (' num2str(tasknum) ') ...']); tasknum = tasknum+1;
            result = hlp_serialize(evaluate(task));
            disp('done with task; opening back link...');
            
            try
                % send off the result
                idx = find(collector==':',1);
                outconn = Socket(collector(1:idx-1), str2num(collector(idx+1:end)));
                disp('connected; now sending...');
                outconn.setTcpNoDelay(1);
                outconn.setSoTimeout(round(1000*opts.timeout_recv));
                out = DataOutputStream(outconn.getOutputStream());
                out.writeInt(taskid);
                out.writeInt(length(result));
                out.writeBytes(result);
                out.flush();
                outconn.close();
                disp('done.');
                disp('waiting for connections...');
            catch
                e = lasterror; %#ok<LERR>
                if isempty(strfind(e.message,'timed out'))
                    disp(['Exception during result forwarding: ' e.message]); end
            end
            
        catch
            conn.close();
            e = lasterror; %#ok<LERR>
            if ~isempty(strfind(e.message,'EOFException'))
                disp(['cancelled.']);
            elseif isempty(strfind(e.message,'timed out'))
                disp(['Exception during task receive: ' e.message]); 
            end
        end
               
    catch
        e = lasterror; %#ok<LERR>
        if isempty(strfind(e.message,'timed out'))
            disp(['Exception during accept: ' e.message]); end
    end
end

function result = evaluate(task)
% evaluate a task
try
    ans = []; %#ok<NOANS>
    eval([task ';']);
    result = ans; %#ok<NOANS>
catch
    result = lasterror; %#ok<LERR>
end
