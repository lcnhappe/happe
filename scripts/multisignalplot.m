function [fh,lh] = multisignalplot(d,varargin)
% -------- [fh,lh] = multisignalplot(d,varargin) ---------
% 
% Plots all of the signals in a data matrix in one figure, separated on the y-axis by a
% the mean(max(data)), or by a user-defined value.
%
%               >>> INPUTS >>>
% Required:
%   d = data matrix
% Optional:
%   Fs = sampling rate...if provided, will plot relative to time on x-axis
%   format = 'r','c'...lets the function know if your data is in row or
%       column format...by default assumes 'c'.
%   col = line color (default 'k');
%   maxval = value to separate lines by...default is mean(max(d)); 
%
%               <<< OUTPUTS <<<
%   fh = handle to axis
%   lh = handle to each line separately
%
% By JMS, 11/04/2015
% ----------------------------------------------------------------------------

if nargin>1; Fs = varargin{1}; 
else Fs = []; end
if nargin>2 && ~isempty(varargin{2}); format = varargin{2};
else format = 'c'; end
if nargin>3 && ~isempty(varargin{3}); col = varargin{3};
else col = 'k'; end
if nargin>4 && ~isempty(varargin{4}); maxval = varargin{4};
else maxval = []; end

% convert to column if row format
if strcmp(format,'r')
    d = d';
end

% get time vec if supplied sampling rate
if ~isempty(Fs)
    time = (1:size(d,1))/Fs;
else
    time = (1:size(d,1));
end

n = size(d,2); 
if isempty(maxval)
    maxval = nanmean(max(d)); % mean of maximum value
end

hold on;
for i = 1:n
    handle(i) = plot(time,d(:,i)-maxval*(i-1),col);
end
fh = gca;
set(fh,'box','off','tickdir','out');
axis tight
ylim([-maxval*(n) maxval]);
    
end