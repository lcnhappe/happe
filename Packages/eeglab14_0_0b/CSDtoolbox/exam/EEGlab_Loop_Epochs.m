% EEGlab_Loop_Epochs.m Sample script to demonstrate three different methods
% of applying a CSD transform to a 3-D data matrix used in EEGlab (i.e.,
% a <channels-by-samples-by-epochs> data matrix). All methods yield
% identical results and require about the same computing time.
%
% 15 July 2010 John J.B. Allen <John.JB.Allen@Arizona.edu>
%
% Updated: $Date: 2010/07/20 10:50:00 $ $Author: jk $
%        - added pre-allocation of memory (improve performance)
%        - use single precision to reduce memory demand
%        - added comments and conclusions

%% claim memory to speed computations
data = single(repmat(NaN,size(EEG.data))); % use single data precision

%% Instruction set #1: Looping method of Jenny
tic                                        % stopwatch on
for ne = 1:length(EEG.epoch)               % loop through all epochs
    myEEG = single(EEG.data(:,:,ne));      % reduce data precision to reduce memory demand
    MyResults = CSD(myEEG,G,H);            % compute CSD for <channels-by-samples> 2-D epoch
    data(:,:,ne) = MyResults;              % assign data output
end
looping_CSD_final = double(data);          % final CSD data
looping_time = toc                         % stopwatch off

data(:,:,:) = NaN;                         % re-initialize data output

%% Instruction set #2: Reshape method of John
tic                                        % stopwatch on
data = CSD(reshape(single(EEG.data), ...   % compute CSD for reshaped data (i.e., use a ...
  EEG.nbchan,EEG.trials*EEG.pnts),G,H);    % <channels-by-(samples*trials)> 2-D data matrix)
data = reshape(data,EEG.nbchan, ...        % reshape CSD output and re-assign to EEGlab ...
  EEG.pnts, EEG.trials);                   % <channels-by-samples-by-epochs> data matrix
reshaped_CSD_final = double(data);         % final CSD data
reshape_time = toc                         % stopwatch off

bigfatzeros = looping_CSD_final ...        % check for differences
                - reshaped_CSD_final;
[min(min(min(bigfatzeros))) ...            % should be [0 0] for identical CSD results
 max(max(max(bigfatzeros)))]            

data(:,:,:) = NaN;                         % re-initialize data output

%% Instruction set #3: TEST 3d matrix as input (should bomb, but works fine)
tic                                        % stopwatch on
data = CSD(single(EEG.data),G,H);          % compute CSD for <channels-by-samples-by-epochs> 3-D data matrix
CSD_3D_final = double(data);               % final CSD data
CSD_3D_time = toc                          % stopwatch off

bigfatzeros = CSD_3D_final ...             % check for differences
                - reshaped_CSD_final;
[min(min(min(bigfatzeros))) ...            % should be [0 0] for identical CSD results
 max(max(max(bigfatzeros)))]            

%% Conclusions
%
% 1) All three methods produce identical results in virtually the same time
%    (note that pre-allocating memory for the output <data> matrix removes
%    the need to sequentially update [i.e., rebuild] a new <data> matrix
%    within a loop statement).
%
% 2) The third method (3-D matrix) has the least convoluted syntax. The
%    reason it works in the first place is that Matlab realizes a 3-D matrix
%    as a sequence of consecutive 2-D matrices (see answer to question #14 
%    on the CSD toolbox FAQ page at URL
%    "http://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/faq.html").
