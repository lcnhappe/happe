function list_properties = single_epoch_channel_properties(EEG,epoch_num,eeg_chans)
if ~isstruct(EEG)
	newdata=EEG;
	clear EEG;
	EEG.data=newdata;
	clear newdata;
end

measure = 1;
% TEMPORAL PROPERTIES

% 1 Median diff value
list_properties(:,measure) = median(diff(EEG.data(eeg_chans,:,epoch_num),[],2),2);
measure = measure + 1;

% 2 Variance of the channels
list_properties(:,measure) = var(EEG.data(eeg_chans,:,epoch_num),[],2);
list_properties(isnan(list_properties(:,measure)),measure)=0;
measure = measure + 1;

% 3 Max difference of each channel
list_properties(:,measure)=(max(EEG.data(eeg_chans,:,epoch_num),[],2)-min(EEG.data(eeg_chans,:,epoch_num),[],2));
measure = measure + 1;

% 4 Deviation from channel mean
list_properties(:,measure)=abs(mean(EEG.data(eeg_chans,:,epoch_num),2)-mean(EEG.data(eeg_chans,:),2));
measure = measure + 1;

for u = 1:size(list_properties,2)
	list_properties(:,u) = list_properties(:,u) - median(list_properties(:,u));
end