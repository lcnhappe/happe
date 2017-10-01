function [results,errors] = hlp_schedule(tasks,varargin)
% Schedule the given tasks across a pool of (possibly remote) workers.
% Results = hlp_schedule(Tasks, Options...)
%
% In:
%   Tasks : cell array of tasks; formatted as
%           * (evaluatable) string
%           * {function_handle, arg1,arg2,arg3, ...}
%           * struct('head',function_handle, 'parts',{{arg1,arg2,arg3, ...}})
%           (see also hlp_beginschedule for further details)
%
%   Options...: optional name-value pairs; see hlp_beginschedule and hlp_endschedule for the options
%
%               'keep': keep this scheduler (and its connections) alive for later re-use (default: true)
%                       if false, the scheduler will be destroyed after use, and re-created during the next run
%
% Out:
%   Results : cell array of results of the scheduled computations (evaluated tasks)
%   Errors  : cell array of exception structs for those results that could not be evaluated (in no particular order)
%
% See also:
%  hlp_worker, hlp_beginschedule, hlp_endschedule
%
% Example:
%  results = hlp_schedule({'sin(randn(10))','exp(randn(10))'},'pool',{'localhost:32547','localhost:32548'})
% 
%                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                       2010-08-29

opts = hlp_varargin2struct(varargin, 'keep',true);

id = hlp_beginschedule(tasks,opts);
[results,errors] = hlp_endschedule(id,opts);
