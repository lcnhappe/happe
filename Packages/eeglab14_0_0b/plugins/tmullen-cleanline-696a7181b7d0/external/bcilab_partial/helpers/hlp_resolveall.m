function outargs = hlp_resolveall(context)
% Look up all symbols in the current hlp_scope scope.
% Scope = hlp_resolveall(Context)
%
% In:
%   Context : Optionally the execution context (stack) at this point, if known.
%             Can be obtained via try/catch.
%
% Out:
%   Scope : a struct with values assigned to symbol names
%
% See also:
%   hlp_resolve
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

% start with the base scope
try
    scope = tracking.stack.base;
catch
    scope = struct();
    tracking.stack.base = scope;
end
% and walk up the stack, overriding existing symbols
for k=frames(end:-1:1)
    frame = tracking.stack.frames.(entries{k}(24:end-14));
    for fn=fieldnames(frame)'
        scope.(fn{1}) = frame.(fn{1}); end
end
outargs = scope;
