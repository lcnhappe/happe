function [results,errors] = hlp_endschedule(sched,varargin)
% Wait for completion of a scheduling operation and return results and errors.
% [Results, Errors] = hlp_endschedule(Id, Options...)
%
% In:
%   Id : scheduler id, obtained from hlp_beginschedule
%
% Out:
%   Results : cell array of results of the scheduled computations (evaluated strings)
%   Errors  : cell array of exception structs for those results that could not be evaluated (in no particular order)
%
% See also:
%   hlp_beginschedule, hlp_worker
%
%                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                       2010-08-29

results = {};
errors = {};

% read options
opts = hlp_varargin2struct(varargin,'keep',false, 'spin_interval',0.1);

if iscell(sched)
    % locally computed results
    for r=1:length(sched)
        results{r} = sched{r}{2}; end
else
    % BLS scheduling
    
    % wait for the scheduler to finish (note: we cannot wait on a condition variable here,
    % as we need the MATLAB thread to be active for managing the reschedule policy)
    while (~sched.done())
        pause(opts.spin_interval); end
    
    % obtain raw results & convert to cell-string array
    strings = sched.results();
    for k=1:length(strings)
        raw{k} = char(strings(k)); end
    
    % terminate scheduler
    if ~opts.keep
        sched.terminate(); end
    
    % evaluate & reorder the string-formatted results
    for r=1:length(raw)
        try
            % evaluate
            [tmp,raw{r}] = evalc(raw{r}); %#ok<ASGLU>
            try
                % put into results
                results{raw{r}{1}} = raw{r}{2};
            catch
                % contains no order id (e.g. was an exception record), put into errors
                errors{end+1} = raw{r};
            end
        catch
            % error evaluating result, put into errors
            errors{end+1} = lasterror; %#ok<LERR>
        end
    end
    
end