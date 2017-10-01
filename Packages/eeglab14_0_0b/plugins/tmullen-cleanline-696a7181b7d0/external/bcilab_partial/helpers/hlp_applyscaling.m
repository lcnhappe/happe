function X = hlp_applyscaling(X, si)
% Apply some previously determined scaling structure
% X = hlp_applyscaling(X, ScaleInfo)
% 
% This is just a convenience tool to implement simple data (e.g. feature) scaling operations.
%
% Examples:
%   scaleinfo = hlp_findscaling(data,'whiten')
%   hlp_applyscaling(data,scaleinfo)
% 
% See also:
%   hlp_findscaling
%
%                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                       2010-03-28

if isfield(si,'add') 
    X = X+repmat(si.add,[size(X,1),1]); end
if isfield(si,'mul') 
    X = X.*repmat(si.mul,[size(X,1),1]); end
if isfield(si,'project') 
    X = X*si.project; end
if isfield(si,{'add','mul','project'}) 
    X(~isfinite(X(:))) = 0; end
