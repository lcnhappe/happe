% eegplugin_cleanline() - EEGLAB plugin for removing line noise
%
% Usage:
%   >> eegplugin_cleanline(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Notes:
%   This plugins consist of the following Matlab files:
%
% Create a plugin:
%   For more information on how to create an EEGLAB plugin see the
%   help message of eegplugin_besa() or visit http://www.sccn.ucsd.edu/eeglab/contrib.html
%
%
% See also: pop_cleanline(), cleanline()

% Copyright (C) 2011 Tim Mullen, SCCN/INC/UCSD
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function vers = eegplugin_cleanline(fig, trystrs, catchstrs)

    vers = 'cleanline';
    if nargin < 3
        error('eegplugin_cleanline requires 3 arguments');
    end;
    
    % add folder to path
    % ------------------
    if exist('cleanline', 'file')
        p = which('eegplugin_cleanline.m');
        p = p(1:findstr(p,'eegplugin_cleanline.m')-1);
        addpath(genpath(p));
    end;
    
    % find import data menu
    % ---------------------
    menu = findobj(fig, 'tag', 'tools');
    
    % menu callbacks
    % --------------
    comcnt = [ trystrs.no_check '[EEG LASTCOM] = pop_cleanline(EEG);' catchstrs.new_and_hist ];
    
    % create menus
    % ------------
    uimenu( menu, 'label', 'CleanLine', 'callback', comcnt,'separator', 'on', 'position',length(get(menu,'children'))+1);
