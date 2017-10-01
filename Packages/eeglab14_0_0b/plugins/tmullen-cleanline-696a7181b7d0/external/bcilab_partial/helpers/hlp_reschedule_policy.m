function result = hlp_reschedule_policy(batchid,inflight,waiting,times,msgs)
% Default re-scheduling policy for the Java-based Scheduler.
% Reschedule = hlp_reschedule_policy(Batch-Id,Inflight-Ids,Event-Times,Event-Messages)
% 
% This function is periodically invoked by the Scheduler with a list of recent events (times in ms and content encoded as strings) and tags of in-flight tasks.
% It is expected to issue the re-scheduling of starved or lost tasks (out of those that are in-flight), depending on some assumptions.
%
% In:
%   Batch-Id     : Identifier of the current batch of tasks (this policy may be invoked for multiple possibly overlapping schedules).
%
%   Inflight-Ids : Cell array of Ids/Tags of in-flight tasks (same as in the Event-Msgs)
%
%   Waiting-Ids : Cell array of Ids/Tags of waiting tasks (same as in the Event-Msgs)
%
%   Event-Times  : Cell array of event timestamps in miliseconds (since beginning of the Scheduler's current batch of tasks).
%
%   Event-Messages: Cell array of event content, indexed like Event-Times.
%
% Out:
%   Reschedule : Java Vector of Integers, referring to the Inflight-Ids that shall be rescheduled (or the empty Vector).
%
%                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                       2010-08-29

% here we keep track of our meta-data for each of the batches
persistent batches;
default_batch = struct('rescheduled',{[]});  % set of already rescheduled tasks for this batch
if isempty(batches)    
    batches = default_batch; % first batch: initialize
elseif length(batches) < batchid    
    batches(batchid) = default_batch; % new batch: append
end

% turn times into seconds relative to the first event
times = [times{:}]/1000;
starttime = min(times);
times = (times-starttime);

% turn inflight into a vector
inflight = [inflight{:}];



% --- begin policies ---
reschedule = [];


% Policy 1: endgame mode: if less than N tasks in flight, re-schedule them to other machines (so that they are being worked on by 2 machines)
if isempty(waiting) && length(inflight) <= 3
    % but do not reschedule them twice
    reschedule = setdiff(inflight,batches(batchid).rescheduled);
end


% Policy 2: identify stragglers, and re-schedule them, assuming that task completion time is normally distributed
try    
    % get a table of recorded tasks
    tasks = struct();
    for m=1:length(msgs)
        msg = msgs{m};
        assigned = isequal(strmatch('assigned',msg),1);
        received = isequal(strmatch('received',msg),1);
        if assigned || received
            tmp = explode(msg,':');
            taskid = ['x' tmp{2}];
            if assigned
                tasks.(taskid).assigned = times(m);
            elseif received
                tasks.(taskid).received = times(m);
            end
            tasks.(taskid).worker = tmp{3};
        end
    end
    
    % get the task completion times for the completed tasks
    completion_times = [];
    for f=fieldnames(tasks)'
        t = tasks.(f{1});
        if isfield(t,{'assigned','received'})
            completion_times(end+1) = t.received - t.assigned; end
    end
    
    % if estimates are reasonable at this point (half of the jobs have been scheduled)
    if length(completion_times) >= length(inflight) && length(completion_times) >= 5
        now = double(tic)/1000000 - starttime;
        % estimate the parameters mu,sigma assuming a truncated normal distribution of completion times
        custompdf = @(x,mu,sigma) normpdf(x,mu,sigma)./normcdf(now,mu,sigma);
        [estim,conf] = mle(completion_times,'pdf',custompdf,'start',[mean(completion_times),std(completion_times)],'lower',[-Inf 0]);      
        mu = estim(1);
        sigma = estim(2);        
        
        % for all inflight (and registered) tasks...
        for i=1:length(inflight)
            t = inflight(i);
            taskid = ['x' num2str(t)];
            % estimate their duration (tasks for which we have no record are assumed to have started before our records began)
            if isfield(tasks,taskid)
                duration(i) = now - tasks.(taskid).assigned;
            else
                duration(i) = now;
            end
        end
        
        % decide which ones we reschedule, based on how long other tasks took so far
        for i=1:length(inflight)
            taskid = inflight(i);
            % duration relative to mean, in std. devs (mahalanobis metric)
            mahal = (duration(i) - mu) / sigma;
            if mahal > 3 && isempty(find(batches(batchid).rescheduled,taskid))
                % is at least >3 std devs, and the task has not yet been rescheduled
                batches(batchid).rescheduled(end+1) = taskid;
                reschedule(end+1) = taskid;
            end
        end
    end
    
catch
    env_handleerror(lasterror);
end


% --- end policies ---

% if ~isempty(reschedule)
%     disp(['Rescheduling job ' num2str(reschedule)]); end

% emit the result as a Java Vector (pre-Generics)
import java.lang.*;
import java.util.*;
result = Vector();
for k=reschedule
    result.add(Integer(k)); end
