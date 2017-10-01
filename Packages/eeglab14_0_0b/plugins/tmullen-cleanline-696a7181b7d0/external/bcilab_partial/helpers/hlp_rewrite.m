function data = hlp_rewrite(data,varargin)
% Rewrite (replace) an input data structure into an output data structure.
% Output = hlp_rewrite(Input,Rules...)
%
% In:
%   Input : some data structure to be rewritten
%
%   Rules... : a comma-separated list of rewrite rules (oldval,newval,oldval,newval,oldval,newval, ...)
%
% Out:
%   Output : the input data structure, but rewritten according to the rules (where they matched).
%
%
% Notes:
%   * No two oldval's may be equal. 
%   * If there is no match, data remains unchanged.
%
% Examples:
%   % rewrite short forms of some string into corresponding long forms
%   hlp_rewrite(myinput, 'hp','highpass', 'lp','lowpass', 'bp','bandpass')
% 
%   % rewrite true to 'on' and false to 'off'
%   hlp_rewrite(myinput, true,'on', false,'off')
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-09-26

old = varargin(1:2:end);
new = varargin(2:2:end);

if iscellstr(old)
    % mapping from strings
    match = strcmp(data,old);
    if any(match)
        data = new{match}; end
else
    % mapping from general structures
    match = cellfun(@(x)isequalwithequalnans(x,data),old);
    if any(match)
        data = new{match}; end        
end
