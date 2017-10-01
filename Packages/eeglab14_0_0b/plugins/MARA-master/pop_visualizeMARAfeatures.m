%   pop_visualizeMARAfeatures() - Display features that MARA's decision 
%       for artifact rejection  is based on 
%                    
%   Usage:
%         >> pop_visualizeMARAfeatures(gcompreject, MARAinfo);
%  
%   Inputs:
%     gcompreject   - array <1 x nIC> containing 1 if component was rejected 
%     MARAinfo      -  struct containing more information about MARA classification 
%                       (output of function <MARA>)
%                           .posterior_artefactprob : posterior probability for each 
%                               IC of being an artefact according to 
%                           .normfeats : <6 x nIC > features computed by MARA for each IC, 
%                                normalized by the training data 
%                           The features are: (1) Current Density Norm, (2) Range
%                           in Pattern, (3) Local Skewness of the Time Series, 
%                           (4) Lambda, (5) 8-13 Hz, (6) FitError. 
%  
% See also: MARA(), processMARA(), pop_selectcomps_MARA()

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
function pop_visualizeMARAfeatures(gcompreject, MARAinfo)

%try
    set(0,'units','pixels');
    resolution = get(0, 'Screensize');
    width = resolution(3);
    height = resolution(4);
    panelsettings.rowsize = 200;
    panelsettings.columnsize = 200;
    panelsettings.columns = floor(width/(2*panelsettings.columnsize));
    panelsettings.rows = floor(height/panelsettings.rowsize);
    panelsettings.numberPerPage = panelsettings.columns * panelsettings.rows;
    panelsettings.pages = ceil(length(gcompreject)/ panelsettings.numberPerPage);

    % display components on a number of different pages 
    for page=1:panelsettings.pages
        selectcomps_1page(page, panelsettings, gcompreject, MARAinfo);    
    end;
   
%catch
%    eeglab_error
%end

 
function EEG = selectcomps_1page(page, panelsettings, gcompreject, MARAinfo)

try, icadefs; 
catch, 
	BACKCOLOR = [0.8 0.8 0.8];
	GUIBUTTONCOLOR   = [0.8 0.8 0.8]; 
end;

% set up the figure
% %%%%%%%%%%%%%%%%
if ~exist('fig')
    mainFig = figure('name', 'Visualize MARA features', 'numbertitle', 'off');
    set(mainFig, 'Color', BACKCOLOR)

    set(gcf,'MenuBar', 'none');
    pos = get(gcf,'Position');
    set(gcf,'Position', [20 20 panelsettings.columnsize*panelsettings.columns panelsettings.rowsize*panelsettings.rows]);
end;

% compute range of components to display
if page < panelsettings.pages
    range = (1 + (panelsettings.numberPerPage * (page-1))) : (panelsettings.numberPerPage + ( panelsettings.numberPerPage * (page-1))); 
else
    range = (1 + (panelsettings.numberPerPage * (page-1))) : length(gcompreject);
end

% draw each component
% %%%%%%%%%%%%%%%%
for i = 1:length(range)
    subplot(panelsettings.rows, panelsettings.columns, i) 
    for j = 1:6
        h = barh(j, MARAinfo.normfeats(j, range(i))); 
        hold on;
        if j <= 4 && MARAinfo.normfeats(j, range(i)) > 0 
            set(h, 'FaceColor', [0.4 0 0]); 
        end
        if j > 4 && MARAinfo.normfeats(j, range(i)) < 0 
            set(h, 'FaceColor', [0.4 0 0]); 
        end        
    end
    axis square;
    if mod(i, panelsettings.columns) == 1
        set(gca,'YTick', 1:6, 'YTickLabel', {'Current Density Norm', ...
            'Range in Pattern', 'Local Skewness', 'lambda', '8-13Hz', 'FitError'})
    else 
        set(gca,'YTick', 1:6, 'YTickLabel', cell(1,6))
    end
    if gcompreject(range(i)) == 1
        title(sprintf('IC %d, p-artifact = %1.2f', range(i),MARAinfo.posterior_artefactprob(range(i))),...
            'Color', [0.4 0 0]); 
        set(gca, 'Color', [1 0.7 0.7])
        %keyboard
    else
        title(sprintf('IC %d, p-artifact = %1.2f', range(i),MARAinfo.posterior_artefactprob(range(i))));
    end
end
  


