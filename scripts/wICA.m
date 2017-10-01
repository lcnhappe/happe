%function [wIC,A,W,IC] = wICA(data, type= 'fastica', plotting= 1, Fs= 250, L=5)
function [wIC,A,W,IC] = wICA(EEG,varargin)
%--------------- function [wIC,A,W] = wICA(data,varargin) -----------------
%
% Performs ICA on data matrix (row vector) and subsequent wavelet
% thresholding to remove low-amplitude activity from the computed ICs.
% This is useful for extracting artifact-only ICs in EEG (for example), and
% then subtracting the artifact-reconstruction from the original data. 
%
%               >>> INPUTS >>>
% Required: 
%   data = data matrix in row format
% Optional:
%   type = "fastica" or "radical"...two different ICA algorithms based on
%       entropy. "fastica" (default) is parametric, "radical" is nonparametric.
%   mult = threshold multiplier...multiplies the computed threshold from
%       "ddencmp" by this number. Higher thresh multipliers = less
%       "background" (or low amp. signal) is kept in the wICs.
%   plotting = 1 or 0. If 1, plots wIC vs. non-wavelet thresholded ICs
%   Fs = sampling rate, (for plotting...default = 1);
%   L = level set for stationary wavelet transform. Higher levels give
%       better frequency resolution, but less temporal resolution. 
%       Default = 5
%   wavename = wavelet family to use. type "wavenames" to see a list of
%       possible wavelets. (default = "coif5");
%
%               <<< OUTPUTS <<<
%   wIC = wavelet-thresholded ICs
%   A = mixing matrix (inv(W)) (optional)
%   W = demixing matrix (inv(A)) (optional)
%   IC = non-wavelet ICs (optional)
%   
%       * you can reconstruct the artifact-only signals as:
%               artifacts = A*wIC;
%       - upon reconstruction, you can then subtract the artifacts from your
%       original data set to remove artifacts, for instance.
%
% Example:
%  n = rand(10,1000);
%  a = [zeros(1,400),[.5,.8,1,2,2.4,2.5,3.5,5,6.3,6,4,3.2,3,1.7,1,-.6,-2.2,-4,-3.6,-3,-1,0],zeros(1,578)];
%  data = n + linspace(0,2,10)'*a;
%  [wIC,A] = wICA(data,[],5,1);
%  ahat = A*wIC;
%  nhat = data-ahat;
%  err = sum(sqrt((nhat-n).^2));

% By JMS, 11/10/2015
%---------------------------------------------------------------------------------------

% check inputs
if nargin>1 && ~isempty(varargin{1})
type=varargin{1}; else type='runica';end
if nargin>2 && ~isempty(varargin{2})
mult=varargin{2};else mult=1;end
if nargin>3 && ~isempty(varargin{3})
plotting=varargin{3}; else plotting=0;end
if nargin>4 && ~isempty(varargin{4})
Fs=varargin{4};else Fs=1;end
if nargin>5 && ~isempty(varargin{5})
L=varargin{5}; else L=5;end
if nargin>6 && ~isempty(varargin{6})
wavename=varargin{6}; else wavename='coif5';end

% run ICA using "runica" or "radical"
if strcmp(type,'runica')
    [OUTEEG, com] = pop_runica(EEG, 'extended',1,'interupt','on'); %runica for parametric, default extended for finding subgaussian distributions
    W = OUTEEG.icaweights*OUTEEG.icasphere;
    A = inv(W);
    IC=reshape(OUTEEG.icaact, size(OUTEEG.icaact,1), []);
    %com = pop_export(OUTEEG,'ICactivationmatrix','ica','on','elec','off','time','off','precision',4);
    %IC = ICactivationmatrix;
elseif strcmp(type,'radical')
    [IC,W] = radical(data); % radical ICA for non-parametric
    A = inv(W);
end

% padding data for proper wavelet transform...data must be divisible by
% 2^L, where L = level set for the stationary wavelet transform
modulus = mod(size(IC,2),2^L); %2^level (level for wavelet)
if modulus ~=0
    extra = zeros(1,(2^L)-modulus);
else
    extra = [];
end
      
% loop through ICs and perform wavelet thresholding
disp('Performing wavelet thresholding');
for s = 1:size(IC,1)
    if ~isempty(extra)
        sig = [IC(s,:),extra]; % pad with zeros
    else
        sig = IC(s,:);
    end
    [thresh,sorh,~] = ddencmp('den','wv',sig); % get automatic threshold value
    thresh = thresh*mult; % multiply threshold by scalar
    swc = swt(sig,L,wavename); % use stationary wavelet transform (SWT) to wavelet transform the ICs
    Y = wthresh(swc,sorh,thresh); % threshold the wavelet to remove small values
    wIC(s,:) = iswt(Y,wavename); % perform inverse wavelet transform to reconstruct a wavelet IC (wIC)
    clear y sig thresh sorh swc 
end

% remove extra padding
if ~isempty(extra)
    wIC = wIC(:,1:end-numel(extra));
end

% plot the ICs vs. wICs
if plotting>0
    disp('Plotting');
    subplot(3,1,1);
        multisignalplot(IC,Fs,'r');
        title('ICs');
    subplot(3,1,2);
        multisignalplot(wIC,Fs,'r');
        title('wICs')
    subplot(3,1,3);
        multisignalplot(IC-wIC,Fs,'r');
        title('Difference (IC - wIC)');
end

end