function sched = hlp_beginschedule(tasks,varargin)
% Begin the scheduling of some set of tasks across a pool of (possibly remote) workers.
% Id = hlp_beginschedule(Tasks,Options...)
%
% Returns a scheduler handle to wait for and obtain results upon completion.
%
% In:
%   Tasks : cell array of tasks; each cell can be formatted as 
%           * (evaluatable) string
%           * {function_handle, arg1,arg2,arg3, ...}
%           * struct('head',function_handle, 'parts',{{arg1,arg2,arg3, ...}})
%
%   Options...: optional name-value pairs, with possible names:
%               'engine': parallelization engine to use, can be one of:
%                         'global': select the global BCILAB setting (bcilab.parallel.engine)
%                         'ParallelComputingToolbox': use the Mathworks Parallel Computing Toolbox (tm); uses resources allocated via the matlabpool command or a configuration file
%                         'BLS': use the BCILAB Scheduler (uses the resources specified in the pool argument) (default)
%                         'Reference': local reference implementation for testing BLS (using the same task serialization mechanism)
%                         'local': do all computations locally, skipping serialization
%
%               'pool': pool of workers to consider for the BLS scheduler (default: {'localhost:23547','localhost:23548', ..., 'localhost:23554'})
%                       if 'global', the global BCILAB setting (bcilab.parallel.pool) will be chosen
%                       (with the BLS engine, an empty pool implies local computation)
%
%               'policy': scheduling policy function for the BLS (default: 'hlp_reschedule_policy')
%                         if 'global', the global BCILAB setting (bcilab.parallel.policy) will be chosen
%
% Out:
%   Id : scheduling id; used to receive results
%
% See also:
%   hlp_endschedule, hlp_worker
%
% Example:
%   id = hlp_beginschedule({'sin(randn(10))','exp(randn(10))'}, 'pool',{'192.168.1.1:23547','192.168.1.2:23547','192.168.1.2:23548','192.168.1.3:23547'});
%   ... optionally do something in the meantime
%   results = hlp_endschedule(id);
%
% Expert note:
%  The 'keep' option of hlp_schedule is also available here and in hlp_endschedule, but it is much harder to use correctly:
%   * if passed as true to hlp_beginschedule, it *must* also be passed as true to hlp_endschedule
%   * nested schedules are not allowed if they use the same worker pool
%
%                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                       2010-08-29

% read options
opts = hlp_varargin2struct(varargin, ...
     'engine','BLS', ...
     'pool', {'localhost:23547','localhost:23548','localhost:23549','localhost:23550','localhost:23551','localhost:23552','localhost:23553','localhost:23554'}, ...
     'policy', 'hlp_reschedule_policy', ...
     'receiver_backlog', 5, ...
     'receiver_timeout', 1, ...
     'reschedule_interval',10, ...
     'keep', false ...
    );

% retrieve global settings, if requested
global bcilab;
if strcmp(opts.engine,'global')
    opts.engine = bcilab.parallel.engine; end
if strcmp(opts.pool,'global')
    opts.pool = bcilab.parallel.pool; end
if strcmp(opts.policy,'global')
    opts.policy = bcilab.parallel.policy; end

% fall back to local computation, if necessary
if strcmp(opts.engine,'BLS') && isempty(opts.pool)
    opts.engine = 'local'; end

if strcmp(opts.engine,'BLS')
    % create a scheduler (Java code, see dependencies/Scheduling-*)
    if opts.keep
        sched = hlp_microcache('schedulers',@(varargin)Scheduler(varargin{:}),opts.pool,opts.policy,opts.receiver_backlog,round(1000*opts.receiver_timeout),round(1000*opts.reschedule_interval));
    else
        sched = Scheduler(opts.pool,opts.policy,opts.receiver_backlog,round(1000*opts.receiver_timeout),round(1000*opts.reschedule_interval));
    end
end

% construct tasks (we augment them by their order id, so we can reorder results properly after they have been collected)
if strcmp(opts.engine,'Reference') || strcmp(opts.engine,'BLS')
    for t=1:length(tasks)
        % task given as Mathematica-style expression struct (see expressions/exp_*)
        if isfield(tasks{t},{'head','parts'})
            tasks{t} = [{tasks{t}.head} tasks{t}.parts]; end
        % task given as {function_handle, arg1, arg2, ...}
        if iscell(tasks{t}) && isa(tasks{t}{1},'function_handle')
            args = hlp_serialize(tasks{t}(2:end));
            tasks{t} = ['feval(',hlp_serialize(tasks{t}{1}), ',', args(2:end-1) ')'];
        end
        % task given as a string
        tasks{t} = ['{' num2str(t) ', ' tasks{t} '}'];
    end
end

% submit tasks
switch opts.engine
    case 'Reference'
        for t=1:length(tasks)
            [tmp,sched{t}] = evalc(tasks{t});  %#ok<ASGLU>
        end
    case 'BLS'
        sched.submit(tasks);
    case 'ParallelComputingToolbox'
        % note: if the interpreter gives an error (pre-2007b), comment it out (the BLS scheduler is a sufficient alternative)
        sched = {};
        parfor t=1:length(tasks)
            if isfield(tasks{t},{'head','parts'})
                sched(t) = {{t,tasks{t}.head(tasks{t}.parts{:})}};
            elseif iscell(tasks{t}) && isa(tasks{t}{1},'function_handle')
                sched(t) = {{t,tasks{t}{1}(tasks{t}{2:end})}};
            else
                sched(t) = {t,eval(tasks{t})};  %#ok<PFBFN>
            end
        end
    case 'local'
        for t=1:length(tasks)
            if isfield(tasks{t},{'head','parts'})
                sched{t} = {t,tasks{t}.head(tasks{t}.parts{:})};
            elseif iscell(tasks{t}) && isa(tasks{t}{1},'function_handle')
                sched{t} = {t,tasks{t}{1}(tasks{t}{2:end})};
            else
                sched{t} = {t,eval(tasks{t})};
            end
        end
    otherwise
        error('unsupported parallelization engine selected.');
end
