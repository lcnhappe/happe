function result = hlp_iscaller(func, level)
% Test whether some function is calling this function at some level(s) of indirection.
% Result = hlp_iscaller(Function, Levels)
%
% It can be specified what levels of nesting outside the function which runs hlp_iscaller are considered.
%
% In:
%   Function : function handle to a function to be tested
%   Levels   : nesting level(s) that shall be tested (default: all)
%              level 1 is the function that invokes hlp_iscaller
%
% Out: 
%   Result   : whether the current code is called by the Caller
%
% Examples:
%   % in a function, test if the calling function, or the caller of the calling function, is named
%   % 'somefunction'
%   hlp_iscaller('somefunction',1:2)
%
% See also:
%   hlp_getcaller
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-04-14

try
    throw; %#ok<LTARG> % fastest way to get an exception
catch context
    if ~exist('level','var')
        level = 2:length(context.stack); end
    level = level+1;
    level(level < 1 | level > length(context.stack)) = [];
    result = any(strcmp(char(func),{context.stack(level).name}));
end
