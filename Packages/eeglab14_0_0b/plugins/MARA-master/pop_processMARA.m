% pop_processMARA() - graphical interface to select MARA's actions
%
% Usage:
%   >> [ALLEEG,EEG,CURRENTSET,com] = pop_processMARA(ALLEEG,EEG,CURRENTSET );
%
% Inputs and Outputs: 
%   ALLEEG      - array of EEG dataset structures
%    EEG        - current dataset structure or structure array
%                  (EEG.reject.gcompreject will be updated)  
%    CURRENTSET - index(s) of the current EEG dataset(s) in ALLEEG
%
% 
% Output:
%   com  - last command, call to itself
%
% See also: processMARA(), pop_selectcomps_MARA(), MARA()

% Copyright (C) 2013 Irene Winkler and Eric Waldburger
% Berlin Institute of Technology, Germany 
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


function [ALLEEG,EEG,CURRENTSET,com] = pop_processMARA (ALLEEG,EEG,CURRENTSET )
com = 'pop_processMARA ( ALLEEG,EEG,CURRENTSET )';

try, icadefs; 
catch, 
	BACKCOLOR = [0.8 0.8 0.8];
	GUIBUTTONCOLOR   = [0.8 0.8 0.8]; 
end;

% set up the figure
% %%%%%%%%%%%%%%%%%
figure('name', 'MARA', ...
    'numbertitle', 'off', 'tag', 'ADEC - Init');
set(gcf,'MenuBar', 'none', 'Color', BACKCOLOR);
pos = get(gcf,'Position');
set(gcf,'Position', [pos(1) 350 350 350]);

if ~strcmp(get(gcf, 'tag'), 'ADEC - Init');
    disp('Aborting plot');
    return;
end;
% set standards for options and store them with guidata()
% order: filter, run_ica, plot data, remove automatically
options = [0, 0, 0, 0, 0];
guidata(gcf, options);

% text
% %%%%%%%%%%%%%%%%%
text = uicontrol(gcf, 'Style', 'text', 'Units','Normalized', 'Position',...
    [0 0.85 0.75 0.06],'BackGroundColor', BACKCOLOR, 'FontWeight', 'bold', ...
    'String', 'Select preprocessing operations:');

text = uicontrol(gcf, 'Style', 'text', 'Units','Normalized', 'Position',...
    [0 0.55 0.75 0.06],'BackGroundColor', BACKCOLOR, 'FontWeight', 'bold', ...
    'String', 'After MARA classification:');                       

% checkboxes
% %%%%%%%%%%%%%%%%%
filterBox = uicontrol(gcf, 'Style', 'checkbox', 'Units','Normalized', 'Position',...
    [0.075 0.75 0.75 0.075], 'String', 'Filter the data', 'Callback', ... 
    'options = guidata(gcf); options(1) = mod(options(1) + 1, 2); guidata(gcf,options);');

icaBox = uicontrol(gcf, 'Style', 'checkbox', 'Units','Normalized', 'Position',...
    [0.075 0.65 0.75 0.075], 'String', 'Run ICA', 'Callback', ...
    'options = guidata(gcf); options(2) = mod(options(2) + 1, 2); guidata(gcf,options);');

% radio button group
% %%%%%%%%%%%%%%%%%
vizBox = uicontrol(gcf, 'Style', 'checkbox', 'Units','Normalized', 'Position',...
    [0.2 0.2 0.9 0.2], 'String', 'Visualize Classifcation Features', 'Enable', 'Off', ...
    'tag', 'vizBox');

h = uibuttongroup('visible','off','Units','Normalized','Position',[0.075 0.15 0.9 0.35]);

% Create two radio buttons in the button group.
radNothing = uicontrol('Style','radiobutton','String','Continue using EEGLab functions',...
    'Units','Normalized','Position',[0.02 0.8 0.9 0.2],'parent',h,'HandleVisibility','off','tag', 'radNothing');
radPlot = uicontrol('Style','radiobutton','String','Plot and select components for removal',...
    'Units','Normalized','Position',[0.02 0.5 0.9 0.2],'parent',h,'HandleVisibility','off',...
    'tag', 'radPlot'); 
radAuto = uicontrol('Style','radiobutton','String','Automatically remove components',...
    'Units','Normalized','Position',[0.02 0.1 0.9 0.2],'parent',h,'HandleVisibility','off','tag', 'radAuto');

set(h,'SelectedObject',radNothing,'Visible','on', 'tag', 'h','BackGroundColor', BACKCOLOR);
set(h, 'SelectionChangeFcn',['s = guihandles(gcf); if get(s.h, ''SelectedObject'') == s.radPlot; ' ...
    'set(s.vizBox,''Enable'', ''on''); else; set(s.vizBox,''Enable'', ''off''); end;']); 

% bottom buttons
% %%%%%%%%%%%%%%%%%
cancel = uicontrol(gcf, 'Style', 'pushbutton', 'Units','Normalized', 'Position',...
    [0.075 0.05 0.4 0.08],'String', 'Cancel','BackgroundColor', GUIBUTTONCOLOR,...
    'Callback', 'close(gcf);'); 

ok = uicontrol(gcf, 'Style', 'pushbutton', 'Units','Normalized', 'Position',...
    [0.5 0.05 0.4 0.08],'String', 'Ok','BackgroundColor', GUIBUTTONCOLOR, ...
    'Callback', ['options = guidata(gcf);' ...
    's = guihandles(gcf); if get(s.h, ''SelectedObject'') == s.radPlot; options(3) = 1; end; '...
    'options(4) = get(s.vizBox, ''Value''); ' ...
    'if get(s.h, ''SelectedObject'') == s.radAuto; options(5) = 1; end;' ...
    'close(gcf); pause(eps); [ALLEEG, EEG, CURRENTSET] = processMARA(ALLEEG,EEG,CURRENTSET, options);']); 
