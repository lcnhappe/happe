function f = hlp_makefunction(decl,workspace)
% Create a function handle from code and a workspace of variables. For hlp_serialize().
% Function = hlp_makefunction(Declaration,Workspace)
% 
% In: 
%   Declaration : The function's declaration string, e.g. '@(x) x*A + B'
%
%   Workspace   : Struct that holds the (non-parameter) variables being referenced in the declaration
%                 (A and B in the above example)
%
% Out:
%   Function : A function handle that can subsequently be evaluated.
%
% Examples:
%   f = hlp_makefunction('@(x)x*A+B',struct('A',2,'B',10));
%   f(2) --> yields 14
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-08-30

% create workspace
for fn=fieldnames(workspace)'
    eval([fn{1} ' = workspace.(fn{1}) ;']); end

% evaluate declaration
f = eval(decl);