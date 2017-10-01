function varargout = hlp_getresult(idx,f,varargin)
% Returns the Result-Idx's output of the function, given the supplied arguments.
% Results... = hlp_getresult(Result-Idx, Function, Arguments...)
%
% In:
%   Result-Idx : index of the result (output) of the given function application; can also 
%                be a vector of indices, then re-emitting those outputs as outputs.
%                note: if this is passed as a cell array, the outputs will be wrapped into a cell 
%                      array.
%
%   Function : function to apply to the Arguments
%
%   Arguments... : list of arguments to the function
%
% Out:
%   Results... : one or more outputs, selected according to Result-Idx from the outputs of 
%                Function(Arguments...), and optionally wrapped into a cell array.
%
% Examples:
%   % get the second output of the sort() function, when applied to some data
%   hlp_getresult(2,@sort,[1 4 2 5 2 1 0 6])
%
% See also:
%   hlp_wrapresults
%
%				          		 Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-28

if isnumeric(idx)
    [tmp{1:max(idx)}] = f(varargin{:});
    varargout = tmp(idx);
elseif iscell(idx)
    [tmp{1:max([idx{:}])}] = f(varargin{:});
    varargout = {tmp([idx{:}])};
else
    error('unsupported index format.');
end