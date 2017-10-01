function list_properties = epoch_properties(EEG,eeg_chans)
list_properties = [];

if length(size(EEG.data)) < 3
	fprintf('Not epoched.\n');
	return;
end

measure = 1;

means = mean(EEG.data(eeg_chans,:),2);

% 1 Epoch's mean deviation from channel means.
for u = 1:size(EEG.data,3)
	list_properties(u,measure) = mean(abs(squeeze(mean(EEG.data(eeg_chans,:,u),2)) - means));
end
measure = measure + 1;

% 2 Epoch variance
list_properties(:,measure) = mean(squeeze(var(EEG.data(eeg_chans,:,:),0,2)));
measure = measure + 1;

% 3 Max amplitude difference
for t = eeg_chans
	for u = 1:size(EEG.data,3)
		ampdiffs(t,u) = max(EEG.data(t,:,u)) - min(EEG.data(t,:,u));
	end
end
list_properties(:,measure) = mean(ampdiffs,1);
measure = measure + 1;

for v = 1:size(list_properties,2)
	list_properties(:,v) = list_properties(:,v) - median(list_properties(:,v));
end