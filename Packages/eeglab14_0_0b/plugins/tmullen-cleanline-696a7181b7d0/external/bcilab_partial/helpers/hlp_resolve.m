function outargs = hlp_resolve(x,default,context)
% Look up a symbol from the current hlp_scope scope.
% Data = hlp_resolve(Symbol,Default,Context)
%
% More precisely, perform a lookup for a symbol in the current (dynamic) scope formed by (possibly
% nested) hlp_scope(s). Returns a cell array of all the values associated with the symbol. Used to
% implement dynamic scoping across functions, used by parts of the expression & online system.
%
% In:
%   Symbol : a symbol (string or function handle) to look up; if omitted, the entire symbol frame
%            will be returned as a struct
%
%   Default : the default value, if the symbol does not exist (note: you may also omit the default
%             if you prefer to get an error in case the symbol does not exist)
%
%   Context : Optionally the execution context (stack) at this point, if known.
%             Can be obtained via try/catch.
%
% Out:
%   Data : the data associated with the symbol, or the symbol itself, if the lookup failed
%
% Examples:
%   % resolve the value of the symbol 'test' against the current dynamic scope
%   value = hlp_resolve('test')
%
% See also:
%   hlp_scope
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-05-03

global tracking;

try
    % obtain the current stack frame entries from the context
    entries = {context.stack.name};
catch context
    % or get the entries from the exception that was just thrown
    entries = {context.stack.name};
end

% check for all entries whether they refer to a stack frame...
marker = 'make_func/@(f,a,frame__';
frames = find(strncmp(entries,marker,23));

% the symbol of interest...
symbol = char(x);
for k=frames
    frameid = entries{k}(24:end-14);
    % ...to look up the symbol in question
    if isfield(tracking.stack.frames.(frameid),symbol)
        outargs = tracking.stack.frames.(frameid).(symbol);
        return;
    end
end

try
    if isfield(tracking.stack.base,symbol)
        % symbol is fetched from the base workspace
        outargs = tracking.stack.base.(symbol);
    else
        % symbol is not in the base workspace: remains unsubstituted
        outargs = default;
    end
catch
    % base workspace doesn't exist yet: create it...
    tracking.stack.base = struct();
    outargs = default;
end
