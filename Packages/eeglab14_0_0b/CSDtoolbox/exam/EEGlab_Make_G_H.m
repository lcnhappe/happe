% EEGlab_Make_G_H.m Sample script to demonstrate how to make G and H matrices with
% CSDToolbox for users of EEGlab.  This specific script loads a neuroscan CNT file
% but this can be edited for other data import options (e.g. EDF or BDF).
% The script opens a typical data file and strips non-EEG channels (which 
% will need to be edited for each user) and then creates the 
% montage and creates G and H matrices. 
%
% 15 July 2010 John J.B. Allen <John.JB.Allen@Arizona.edu>
%
% Updated: $Date: 2010/07/19 16:11:00 $ $Author: jk $
%        - added disclaimer that some code may depend on a specific release of EEGlab
%        - revised method of saving G and H matrices for later import
%        - simplified directory references where possible
%        - added some comments

%% Paths and EEGLab
addpath(genpath('C:\Program Files\MATLAB\eeglab8_0_3_5b')); % Point to your version of EEGLab
eeglab;                                                     % Start EEGlab  
addpath(genpath(('C:\Program Files\MATLAB\CSDtoolbox')));   % Point to place where you put CSDToolbox

%% OPEN REPRESENTATIVE DATA FILE
% This instantiation is for neuroscan, but you can also open other formats
% (see sample at end of this script)
[FileName,PathName,FilterIndex] = uigetfile('G:\PhysioData\*.cnt','Choose File to Open');
FullName = [PathName FileName];
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadcnt(FullName, 'dataformat', 'int32', 'keystroke', 'on');  % change to int16 for 16 bit files (older neuroscan)
eeglab redraw;

%% Extract only 64 relevant channels (strip out EKG and EOG)
% Edit to match your montage
EEG = pop_select( EEG,'channel',{'Fp1' 'AF7' 'AF3' 'F1' 'F3' 'F5' 'F7' 'FT7' 'FC5' 'FC3' 'FC1' 'C1' 'C3' 'C5' 'T7' 'TP7' 'CP5' 'CP3' 'CP1' 'P1' 'P3' 'P5' 'P7' 'P9' 'PO7' 'PO3' 'O1' 'Iz' 'Oz' 'POz' 'Pz' 'CPz' 'Fpz' 'Fp2' 'AF8' 'AF4' 'AFz' 'Fz' 'F2' 'F4' 'F6' 'F8' 'FT8' 'FC6' 'FC4' 'FC2' 'FCz' 'Cz' 'C2' 'C4' 'C6' 'T8' 'TP8' 'CP6' 'CP4' 'CP2' 'P2' 'P4' 'P6' 'P8' 'P10' 'PO8' 'PO4' 'O2'});
% Note 1: This code caused an error when running a prior release of EEGlab (version 6.0)
% Note 2: Alternatively, one could directly exclude the unwanted channels by unsing the 'nochannel' argument as in
%         EEG = pop_select( EEG,'nochannel',{'EKG' 'EOG'});
%         which will implicitly keep all other (EEG) channels.
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','BDF file 64 chans','gui','off'); 

%% Get usable list of electrodes from EEGlab data structure
for site = 1:64 
    trodes{site}=(EEG.chanlocs(site).labels);
end;
trodes=trodes';

%% Get Montage for use with CSD Toolbox
Montage_64=ExtractMontage('10-5-System_Mastoids_EGI129.csd',trodes);
MapMontage(Montage_64);

%% Derive G and H!
[G,H] = GetGH(Montage_64);

%% Save G and H to later import when doing the CSD transform on files
% save('G:\PhysioData\MN_Fear\G.mat', 'G');
% save('G:\PhysioData\MN_Fear\H.mat', 'H');

% revised method to store G and H matrices with CSD montage for later import
Montage = Montage_64;                             % use generic variable name
save G:\PhysioData\CSDmontage_64.mat G H Montage; % save variables to Matlab file
clear G H Montage;                                % remove variables from workspace
load G:\PhysioData\CSDmontage_64.mat;             % restore variables to the workspace

%% SAMPLE TO OPEN BIOSEMI -- can use instead of neuroscan version above
% [FileName,PathName,FilterIndex] = uigetfile('L:\Physiodata\*.bdf','Choose File to Open')
% FullName = [PathName FileName];
% [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
% EEG = pop_biosig(FullName, 'ref',[65 66] ,'blockepoch','off'); % Choose ref sites per your montage
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
% eeglab redraw;
