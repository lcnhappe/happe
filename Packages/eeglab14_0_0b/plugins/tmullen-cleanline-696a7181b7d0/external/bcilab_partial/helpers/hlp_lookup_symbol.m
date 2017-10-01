function outargs = hlp_lookup_symbol(x,context)
% Internal. Look up a symbol from the current ult_block / workspace scope.
% Data = hlp_lookup_symbol(Symbol,Context)
%
% More precisely, perform a lookup for a symbol in the current (dynamic) scope formed by (possibly nested) utl_block(s)
% returns a cell array of all the values associated with the symbol.
% Used to implement dynamic scoping across functions, used by parts of the expression system.
%
% In:
%   Symbol : a symbol (function handle or string).
%
%   Context : Optionally the execution context (stack) at this point, if known.
%             Can be obtained via try/catch.
%
% Out:
%   Data : the data associated with the symbol, or the symbol itself, if the lookup failed.
%
%                                        Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                        2010-05-03
global bcilab;

% the symbol of interest...
symbol = char(x);

try
    % obtain the current stack frame entries from context
    entries = {context.stack.name};
catch
    % context did not exist: obtain it
    context = lasterror;  %#ok<LERR>
    % and get entries, too
    entries = {context.stack.name};
end

% check for all entries whether they refer to a stack frame...
for k=find(strncmp(entries,'utl_block/@(x,active_stackframe__',33))
    frameid = entries{k}(34:end-12);
    % ...to look up the symbol in question    
    if isfield(bcilab.stackframes.ov.(frameid),symbol)
        outargs = {bcilab.stackframes.ov.(frameid).(symbol)};
        return;
    end
end

if ~isfield(bcilab.workspace.ov,symbol)
    % symbol is not in the global workspace: remains unsubstituted
    outargs = {x};
else
    % symbol is fetched from the global workspace
    outargs = {bcilab.workspace.ov.(symbol)};
end
