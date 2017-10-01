% eegplugin_MARA() - EEGLab plugin to classify artifactual ICs based on 
%                   6 features from the time domain, the frequency domain, 
%                   and the pattern
%
% Inputs:
%   fig           - [integer]  EEGLAB figure
%   try_strings   - [struct] "try" strings for menu callbacks.
%   catch_strings - [struct] "catch" strings for menu callbacks.
%
% See also:  pop_processMARA(), processMARA(), MARA()

% Copyright (C) 2013  Irene Winkler and Eric Waldburger
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




function eegplugin_MARA( fig, try_strings, catch_strings)

toolsmenu = findobj(fig, 'tag', 'tools');
h = uimenu(toolsmenu, 'label', 'IC Artifact Classification (MARA)');
       


uimenu(h, 'label', 'MARA Classification', 'callback', ...
           [try_strings.no_check  ...
           '[ALLEEG,EEG,CURRENTSET,LASTCOM]= pop_processMARA( ALLEEG ,EEG ,CURRENTSET );' ...
           catch_strings.add_to_hist ]); 
       
uimenu(h, 'label', 'Visualize Components', 'tag', 'MARAviz', 'Enable', ...
            'off', 'callback', [try_strings.no_check  ...
           'EEG = pop_selectcomps_MARA(EEG);  pop_visualizeMARAfeatures(EEG.reject.gcompreject, EEG.reject.MARAinfo); ' ...
           catch_strings.add_to_hist ]); 
       
uimenu(h, 'label', 'About', 'Separator', 'on', 'Callback', ... 
    ['warndlg2(sprintf([''MARA automatizes the process of hand-labeling independent components for ', ...
    'artifact rejection. It is a supervised machine learning algorithm that learns from ', ...
    'expert ratings of 1290 components. Features were optimized to solve the binary classification problem ', ...
    'reject vs. accept.\n \n', ...
    'If you have questions or suggestions about the toolbox, please contact \n ', ...
    'Irene Winkler, TU Berlin irene.winkler@tu-berlin.de \n \n ', ...
    'Reference: \nI. Winkler, S. Haufe, and M. Tangermann, Automatic classification of artifactual', ...
    'ICA-components for artifact removal in EEG signals, Behavioral and Brain Functions, 7, 2011.''])', ...
    ',''About MARA'');']);

