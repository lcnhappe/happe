%   pop_selectcomps_MARA() - Display components with checkbox to label 
%       them for artifact rejection  
%                    
%   Usage:
%         >> EEG = pop_selectcomps_MARA(EEG, gcompreject_old);
%  
%   Inputs:
%     EEG - Input dataset with rejected components (saved in
%           EEG.reject.gcompreject)
%
%  Optional Input: 
%     gcompreject_old - gcompreject to revert to in case the "Cancel"
%     button is pressed
%  
%   Output:
%     EEG - Output dataset with updated rejected components
%  
% See also: processMARA(), pop_selectcomps()

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
function [EEG, com] = pop_selectcomps_MARA(EEG, varargin)


if isempty(EEG.reject.gcompreject)
    EEG.reject.gcompreject = zeros(size(EEG.icawinv,2)); 
end

if not(isempty(varargin))
    gcompreject_old = varargin{1}; 
    com = [ 'pop_selectcomps_MARA(' inputname(1) ',' inputname(2) ');'];
else
    gcompreject_old = EEG.reject.gcompreject; 
    com = [ 'pop_selectcomps_MARA(' inputname(1) ');'];
end


try
    set(0,'units','pixels');
    resolution = get(0, 'Screensize');
    width = resolution(3);
    height = resolution(4);
    panelsettings = {};
    panelsettings.completeNumber = length(EEG.reject.gcompreject);
    panelsettings.rowsize = 250;
    panelsettings.columnsize = 300;
    panelsettings.column = floor(width/panelsettings.columnsize);
    panelsettings.rows = floor(height/panelsettings.rowsize);
    panelsettings.numberPerPage = panelsettings.column * panelsettings.rows;
    panelsettings.pages = ceil(length(EEG.reject.gcompreject)/ panelsettings.numberPerPage);
    panelsettings.incy = 110;
    panelsettings.incx = 110;

    % display components on a number of different pages 
    for page=1:panelsettings.pages
        EEG = selectcomps_1page(EEG, page, panelsettings, gcompreject_old);    
    end;
    
    
    checks = findall(0, 'style', 'checkbox');
    for i = 1: length(checks)
        set(checks(i), 'Enable', 'on')
    end
catch
    eeglab_error
end

 
function EEG = selectcomps_1page(EEG, page, panelsettings, gcompreject_old)
% Display components with checkbox to label them for artifact rejection

try, icadefs; 
catch, 
	BACKCOLOR = [0.8 0.8 0.8];
	GUIBUTTONCOLOR   = [0.8 0.8 0.8]; 
end;


% compute components (for plotting the spectrum) 
if length(size(EEG.data)) == 3
    s = size(EEG.data); 
    data = reshape(EEG.data, [EEG.nbchan, prod(s(2:3))]); 
    icacomps = EEG.icaweights * EEG.icasphere * data;
else 
    icacomps = EEG.icaweights * EEG.icasphere * EEG.data;
end

% set up the figure
% %%%%%%%%%%%%%%%%
if ~exist('fig')
    mainFig = figure('name', [ 'MARA on dataset: ' EEG.setname ''], 'numbertitle', 'off', 'tag', 'ADEC - Plot');
    set(mainFig, 'Color', BACKCOLOR)

    set(gcf,'MenuBar', 'none');
    pos = get(gcf,'Position');
    set(gcf,'Position', [20 20 panelsettings.columnsize*panelsettings.column panelsettings.rowsize*panelsettings.rows]);

    panelsettings.sizewx = 50/panelsettings.column;
    panelsettings.sizewy = 80/panelsettings.rows;
    pos = get(gca,'position'); % plot relative to current axes
    hh = gca;
    panelsettings.q = [pos(1) pos(2) 0 0];
    panelsettings.s = [pos(3) pos(4) pos(3) pos(4)]./100;
    axis off;
end;

% compute range of components to display
if page < panelsettings.pages
    range = (1 + (panelsettings.numberPerPage * (page-1))) : (panelsettings.numberPerPage + ( panelsettings.numberPerPage * (page-1))); 
else
    range = (1 + (panelsettings.numberPerPage * (page-1))) : panelsettings.completeNumber;
end

data = struct; 
data.gcompreject =  EEG.reject.gcompreject;
data.gcompreject_old = gcompreject_old; 
guidata(gcf, data);

% draw each component
% %%%%%%%%%%%%%%%%
count = 1;
for ri = range

    % compute coordinates
    X = mod(count-1, panelsettings.column)/panelsettings.column * panelsettings.incx-10;  
    Y = (panelsettings.rows-floor((count-1)/panelsettings.column))/panelsettings.rows * panelsettings.incy - panelsettings.sizewy*1.3;  

    % plot the head
    if ~strcmp(get(gcf, 'tag'), 'ADEC - Plot');
        disp('Aborting plot');
        return;
    end;
    ha = axes('Units','Normalized', 'Position',[X Y panelsettings.sizewx*0.85 panelsettings.sizewy*0.85].*panelsettings.s+panelsettings.q);      
    topoplot( EEG.icawinv(:,ri), EEG.chanlocs, 'verbose', 'off', 'style' , 'fill');
    axis square;

    % plot the spectrum
    ha = axes('Units','Normalized', 'Position',[X+1.05*panelsettings.sizewx Y-1 (panelsettings.sizewx*1.15)-1 panelsettings.sizewy-1].*panelsettings.s+panelsettings.q);
    [pxx, freq] = pwelch(icacomps(ri,:), ones(1, EEG.srate), [], EEG.srate, EEG.srate);
    pxx = 10*log10(pxx * EEG.srate/2);

    plot(freq, pxx, 'LineWidth', 2)
    xlim([0 50]); 
    grid on; 
    xlabel('Hz')
    set(gca, 'Xtick', 0:10:50)
    if isfield(EEG.reject, 'MARAinfo')
        title(sprintf('Artifact Probability = %1.2f', ...
            EEG.reject.MARAinfo.posterior_artefactprob(ri))); 
    end
    uicontrol(gcf, 'Style', 'checkbox', 'Units','Normalized', 'Position',...
        [X+panelsettings.sizewx*0.5 Y+panelsettings.sizewy panelsettings.sizewx, ...
        panelsettings.sizewy*0.2].*panelsettings.s+panelsettings.q, ...
        'Enable', 'off','tag', ['check_' num2str(ri)], 'Value', EEG.reject.gcompreject(ri), ...
        'String', ['IC' num2str(ri) ' - Artifact?'], ... 
        'Callback', {@callback_checkbox, ri});
    drawnow;
    count = count +1;
end;

% draw the botton buttons
% %%%%%%%%%%%%%%%%
if ~exist('fig')
    % cancel button
    cancelcommand = ['openplots = findall(0, ''tag'', ''ADEC - Plot'');', ...
            'data = guidata(openplots(1)); close(openplots);', ...
            'EEG.reject.gcompreject = data.gcompreject_old;']; 
    hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'Cancel', ...
        'Units','Normalized', 'BackgroundColor', GUIBUTTONCOLOR, ...
            'Position',[-10 -11  15 panelsettings.sizewy*0.25].*panelsettings.s+panelsettings.q, ...
            'callback', cancelcommand );
    okcommand = ['data = guidata(gcf); EEG.reject.gcompreject = data.gcompreject; ' ...
        '[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); '...
        'eegh(''[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);''); '...
        'close(gcf);' ]; 
    % ok button
    hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'OK', ...
        'Units','Normalized', 'BackgroundColor', GUIBUTTONCOLOR, ...
           'Position',[10 -11  15 panelsettings.sizewy*0.25].*panelsettings.s+panelsettings.q, ...
           'callback', okcommand);    
end;




function callback_checkbox(hObject,eventdata, position) 
openplots = findall(0, 'tag', 'ADEC - Plot');
data = guidata(openplots(1)); 
data.gcompreject(position) = mod(data.gcompreject(position) + 1, 2);
%save changes in every plot
for i = 1: length(openplots)
    guidata(openplots(i), data);
end




