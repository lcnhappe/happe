% MARA() - Automatic classification of multiple artifact components
%          Classies artifactual ICs based on 6 features from the time domain, 
%           the frequency domain, and the pattern
%
% Usage:
%   >> [artcomps, info] = MARA(EEG);
%
% Inputs: 
%   EEG         - input EEG structure
% 
% Outputs:
%   artcomps    - array containing the numbers of the artifactual 
%                 components
%   info        - struct containing more information about MARA classification 
%                   .posterior_artefactprob : posterior probability for each 
%                            IC of being an artefact 
%                   .normfeats : <6 x nIC > features computed by MARA for each IC, 
%                            normalized by the training data 
%                      The features are: (1) Current Density Norm, (2) Range
%                      in Pattern, (3) Local Skewness of the Time Series, 
%                      (4) Lambda, (5) 8-13 Hz, (6) FitError. 
%
%  For more information see: 
%  I. Winkler, S. Haufe, and M. Tangermann, Automatic classification of artifactual ICA-components 
%  for artifact removal in EEG signals, Behavioral and Brain Functions, 7, 2011.
%
% See also: processMARA()

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
function [artcomps, info] = MARA(EEG)
try 
    %%%%%%%%%%%%%%%%%%%%
    %%  Calculate features from the pattern (component map)
    %%%%%%%%%%%%%%%%%%%%
    % extract channel labels
    clab = {};
    for i=1:length(EEG.chanlocs)
         clab{i} = EEG.chanlocs(i).labels;
    end

    % cut to channel labels common with training data 
    load('fv_training_MARA'); %load struct fv_tr    
    [clab_common i_te i_tr ] = intersect(upper(clab), upper(fv_tr.clab));
    clab_common = fv_tr.clab(i_tr); 
    if length(clab_common) == 0
        error(['There were no matching channeldescriptions found.' , ... 
        'MARA needs channel labels of the form Cz, Oz, F3, F4, Fz, etc. Aborting.'])
    end
    patterns = (EEG.icawinv(i_te,:));
    [M100 idx] = get_M100_ADE(clab_common); %needed for Current Density Norm

    disp('MARA is computing features. Please wait');
    %standardize patterns
    patterns = patterns./repmat(std(patterns,0,1),length(patterns(:,1)),1);

    %compute current density norm
    feats(1,:) = log(sqrt(sum((M100*patterns(idx,:)).^2)));
    %compute spatial range
    feats(2,:) = log(max(patterns) - min(patterns));

    %%%%%%%%%%%%%%%%%%%%
    %%  Calculate time and frequency features 
    %%%%%%%%%%%%%%%%%%%%
    %compute time and frequency features (Current Density Norm, Range Within Pattern,
    %Average Local Skewness, Band Power 8 - 13 Hz) 
    feats(3:6,:) = extract_time_freq_features(EEG);
    disp('Features ready');


    %%%%%%%%%%%%%%%%%%%%%%
    %%  Adapt train features to clab 
    %%%%%%%%%%%%%%%%%%%%
     fv_tr.pattern = fv_tr.pattern(i_tr, :);
     fv_tr.pattern = fv_tr.pattern./repmat(std(fv_tr.pattern,0,1),length(fv_tr.pattern(:,1)),1);
     fv_tr.x(2,:) = log(max(fv_tr.pattern) - min(fv_tr.pattern));
     fv_tr.x(1,:) = log(sqrt(sum((M100 * fv_tr.pattern).^2))); 

    %%%%%%%%%%%%%%%%%%%%
    %%  Classification 
    %%%%%%%%%%%%%%%%%%%%
    [C, foo, posterior] = classify(feats',fv_tr.x',fv_tr.labels(1,:));
    artcomps = find(C == 0)';
    info.posterior_artefactprob = posterior(:, 1)'; 
    info.normfeats = (feats - repmat(mean(fv_tr.x, 2), 1, size(feats, 2)))./ ...
                    repmat(std(fv_tr.x,0, 2), 1, size(feats, 2)); 
catch
  eeglab_error; 
  artcomps = []; 
end
end

function features = extract_time_freq_features(EEG)
%                             - 1st row: Average Local Skewness
%                             - 2nd row: lambda
%                             - 3rd row: Band Power 8 - 13 Hz 
%                             - 4rd row: Fit Error
%                           
data = EEG.data;
fs = EEG.srate; %sampling frequency

% transform epoched data into continous data
if length(size(data)) == 3
    s = size(data); 
    data = reshape(data, [EEG.nbchan, prod(s(2:3))]); 
end

%downsample (to 100-200Hz) 
factor = max(floor(EEG.srate/100),1); 
data = data(:, 1:factor:end); 
fs = round(fs/factor); 
 
%compute icaactivation and standardise variance to 1
icacomps = (EEG.icaweights * EEG.icasphere * data)';
icacomps = icacomps./repmat(std(icacomps,0,1),length(icacomps(:,1)),1);
icacomps = icacomps';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate featues  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ic=1:length(icacomps(:,1))  %for each component
    fprintf('.');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Proc Spectrum for Channel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [pxx, freq] = pwelch(icacomps(ic,:), ones(1, fs), [], fs, fs);
    pxx = 10*log10(pxx * fs/2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The average log band power between 8 and 13 Hz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = 0;
    for i = 8:13 
        p = p + pxx(find(freq == i,1));
    end
    Hz8_13 = p / (13-8+1);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % lambda and FitError: deviation of a component's spectrum from
    % a protoptypical 1/frequency curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p1.x = 2; %first point: value at 2 Hz
    p1.y = pxx(find(freq == p1.x,1)); 

    p2.x = 3; %second point: value at 3 Hz
    p2.y = pxx(find(freq == p2.x,1));
    
    %third point: local minimum in the band 5-13 Hz
    p3.y = min(pxx(find(freq == 5,1):find(freq == 13,1)));
    p3.x = freq(find(pxx == p3.y,1));

    %fourth point: min - 1 in band 5-13 Hz
    p4.x = p3.x - 1;
    p4.y = pxx(find(freq == p4.x,1));

    %fifth point: local minimum in the band 33-39 Hz
    p5.y = min(pxx(find(freq == 33,1):find(freq == 39,1)));
    p5.x = freq(find(pxx == p5.y,1));
    
    %sixth point: min + 1 in band 33-39 Hz
    p6.x = p5.x + 1;
    p6.y = pxx(find(freq == p6.x,1));
    
    pX = [p1.x; p2.x; p3.x; p4.x; p5.x; p6.x];
    pY = [p1.y; p2.y; p3.y; p4.y; p5.y; p6.y];
    
    myfun = @(x,xdata)(exp(x(1))./ xdata.^exp(x(2))) - x(3);
    xstart = [4, -2, 54];
    try
        fittedmodel = lsqcurvefit(myfun,xstart,double(pX),double(pY), [], [], optimset('Display', 'off'));
    catch
        try
            % If the optimization toolbox is missing we try with the CurveFit toolbox
            opt = fitoptions('Method','NonlinearLeastSquares','Startpoint',xstart);
            myfun = fittype('exp(x1)./x.^exp(x2) - x3;','options',opt);
            fitobject = fit(double(pX),double(pY),myfun);
            fittedmodel = [fitobject.x1, fitobject.x2, fitobject.x3];
        catch
            % If the CurveFit toolbox is also missing we try with the Statistitcs toolbox
            myfun = @(p,xdata)(exp(p(1))./ xdata.^exp(p(2))) - p(3);
            mdl = NonLinearModel.fit(double(pX),double(pY),myfun,xstart);
            fittedmodel = mdl.Coefficients.Estimate(:)';
        end
    end   
    
    %FitError: mean squared error of the fit to the real spectrum in the band 2-40 Hz.
    ts_8to15 = freq(find(freq == 8) : find(freq == 15));
    fs_8to15 = pxx(find(freq == 8) : find(freq == 15)); 
    fiterror = log(norm(myfun(fittedmodel, ts_8to15)-fs_8to15)^2); 
    
    %lambda: parameter of the fit
    lambda = fittedmodel(2); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Averaged local skewness 15s
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    interval = 15; 
    abs_local_scewness = [];
    for i=1:interval:length(icacomps(ic,:))/fs-interval
        abs_local_scewness = [abs_local_scewness, abs(skewness(icacomps(ic, i * fs:(i+interval) * fs)))];
    end
    
    if isempty(abs_local_scewness)
        error('MARA needs at least 15ms long ICs to compute its features.')
    else
        mean_abs_local_scewness_15 = log(mean(abs_local_scewness));
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Append Features 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    features(:,ic)= [mean_abs_local_scewness_15, lambda, Hz8_13, fiterror]; 
end
disp('.');

end

function [M100, idx_clab_desired] = get_M100_ADE(clab_desired)
% [M100, idx_clab_desired] = get_M100_ADEC(clab_desired)
%
% IN  clab_desired - channel setup for which M100 should be calculated
% OUT M100 
%     idx_clab_desired
% M100 is the matrix such that  feature = norm(M100*ica_pattern(idx_clab_desired), 'fro')
%
% (c) Stefan Haufe

lambda = 100;

load inv_matrix_icbm152; %L (forward matrix 115 x 2124 x 3), clab (channel labels) 

[cl_ ia idx_clab_desired] = intersect(clab, clab_desired);
F = L(ia, :, :); %forward matrix for desired channels labels
[n_channels m foo] = size(F);  %m = 2124, number of dipole locations 
F = reshape(F, n_channels, 3*m);

%H - matrix that centralizes the pattern, i.e. mean(H*pattern) = 0
H = eye(n_channels) -  ones(n_channels, n_channels)./ n_channels; 
%W - inverse of the depth compensation matrix Lambda
W = sloreta_invweights(L);

L = H*F*W;

%We have inv(L'L +lambda eye(size(L'*L))* L' = L'*inv(L*L' + lambda
%eye(size(L*L')), which is easier to calculate as number of dimensions is 
%much smaller

%calulate the inverse of L*L' + lambda * eye(size(L*L')
[U D] = eig(L*L');
d = diag(D);
di = d+lambda;
di = 1./di;
di(d < 1e-10) = 0;
inv1 = U*diag(di)*U';  %inv1 = inv(L*L' + lambda *eye(size(L*L'))

%get M100
M100 = L'*inv1*H;
       
end
    
    
function W = sloreta_invweights(LL)
% inverse sLORETA-based weighting
%
% Synopsis:
%   W = sloreta_invweights(LL);
%   
% Arguments:
%   LL: [M N 3] leadfield tensor
%   
% Returns:
%   W: [3*N 3*N] block-diagonal matrix of weights
%
% Stefan Haufe, 2007, 2008
%
% License
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see http://www.gnu.org/licenses/.

[M N NDUM]=size(LL);
L=reshape(permute(LL, [1 3 2]), M, N*NDUM);

L = L - repmat(mean(L, 1), M, 1);

T = L'*pinv(L*L');

W = spalloc(N*NDUM, N*NDUM, N*NDUM*NDUM);
for ivox = 1:N
  W(NDUM*(ivox-1)+(1:NDUM), NDUM*(ivox-1)+(1:NDUM)) = (T(NDUM*(ivox-1)+(1:NDUM), :)*L(:, NDUM*(ivox-1)+(1:NDUM)))^-.5;
end

ind = [];
for idum = 1:NDUM
  ind = [ind idum:NDUM:N*NDUM];
end
W = W(ind, ind);

end



function [i_te, i_tr] = findconvertedlabels(pos_3d, chanlocs)
% IN  pos_3d  - 3d-positions of training channel labels
%     chanlocs - EEG.chanlocs structure of data to be classified

    %compute spherical coordinates theta and phi for the training channel
    %label
    [theta, phi, r] = cart2sph(pos_3d(1,:),pos_3d(2,:), pos_3d(3,:));
    theta = theta - pi/2; 
    theta(theta < -pi) = theta(theta < -pi) + 2*pi; 
    theta = theta*180/pi; 
    phi = phi * 180/pi;
    theta(find(pos_3d(1,:) == 0 & pos_3d(2,:) == 0)) = 0; %exception for Cz

    
    clab_common = {}; 
    i_te = []; 
    i_tr = []; 
    
    %For each channel in EEG.chanlocs, try to find matching channel in
    %training data
    for chan = 1:length(chanlocs)
        if not(isempty(chanlocs(chan).sph_phi))
            idx = find((theta <= chanlocs(chan).sph_theta + 6) ... 
            & (theta >= chanlocs(chan).sph_theta - 6) ...
            & (phi <= chanlocs(chan).sph_phi + 6) ... 
            & (phi >= chanlocs(chan).sph_phi - 6)); 
            if not(isempty(idx))
                    i_tr = [i_tr, idx(1)];
                    i_te = [i_te, chan]; 
            end
        end
    end
end
    