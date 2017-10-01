function a = hlp_wrapresults(f,varargin)
% Wraps all outputs produced by function f for the given arguments into a cell array.
% Results = hlp_wrapresults(Function, Arguments...)
%
% In:
%   Function     : some function handle to execute
%
%   Arguments... : list of arguments to pass to the function
%
% Out:
%   Results : cell array of all function results
%
%
% Notes: 
%   It is not (currently) possible to efficiently determine the number of out-args for a 
%   varargout function; in this case, at most 10 outputs are supported.
%
% Examples:
%   % wrap both outputs of a particualr sort() call into a cell array
%   results = hlp_wrapresults(@sort,[1 4 2 1 3 6 0])
%
% See also:
%   hlp_getresult
%
%					         	 Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-28

% find out how many arguments f can maximally produce
len = nargout(f);
if len < 0
    % for varargout functions, we assume at most 10
    len = 10; end
% ... then invoke the function
while len >= 1
    try
        % get the appropriate number of results from f()
        [a{1:len}] = f(varargin{:});
        return;
    catch e
        % got an exception, check if it is outarg-related
        if ~any(strcmp(e.identifier,{'MATLAB:TooManyOutputs','MATLAB:maxlhs','MATLAB:unassignedOutputs'}))
            % it isn't: rethrow
            rethrow(e); end
        % but if it was, we need to retry with fewer out-arguments
        len = len-1;
    end
end

% len = 0: f produces no outputs
f(varargin{:});
a = {};
