function [EEGclean, commstr, g, Sorig, Sclean, f, amps, freqs] = pop_cleanline(EEG,varargin)

% Mandatory             Information
% --------------------------------------------------------------------------------------------------
% EEG                   EEGLAB data structure
% --------------------------------------------------------------------------------------------------
%
% Optional              Information
% --------------------------------------------------------------------------------------------------
% Type 'doc cleanline' for additional arguments
%
% See Also: cleanline()

% Author: Tim Mullen, SCCN/INC/UCSD Copyright (C) 2011
% Date:   Nov 20, 2011
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
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


% defaults
EEGclean = EEG;
[commstr, Sorig, Sclean, f, amps, freqs, g] = deal([]);


if nargin<2
    % render GUI
    g = arg_guidialog(@cleanline,'Parameters',{'EEG',EEG},'Title','CleanLine Options','Invoke',false);
    
    if isempty(g)
        return;
    end
else
    g = hlp_varargin2struct(varargin);
end


[EEGclean, Sorig, Sclean, f, amps, freqs, g] = cleanline('EEG',EEG,g);

if ~isempty(Sorig)
    
    % plot the original and cleaned spectra
    eegplot(Sorig,'data2',Sclean, ...
        'title',sprintf('Original and Cleaned %s spectra for selected %s',fastif(g.normSpectrum,'normalized',''),g.sigtype),'srate',length(f)/f(end), ...
        'winlength',f(end),'submean','on'); %,'trialstag',1/length(f));
    
    ax = findobj(gcf,'tag','eegaxis');
    xlabel(ax,'Frequency (Hz)');
    title(ax,sprintf('Original and Cleaned %s spectra for selected %s',fastif(g.normSpectrum,'normalized',''),g.sigtype));
    set(ax,'Yticklabel',[{''}; cellstr(num2str(g.chanlist(end:-1:1)'))]);
    %     set(ax,'Xtick',1:10:length(f));
    %     set(ax,'Xticklabel',f(1:10:end));
    plts = get(ax,'children');
    legend([plts(end) plts(1)],'original','cleaned');
    
end

commstr = sprintf('EEG = pop_cleanline(EEG, %s);', vararg2str(hlp_struct2varargin(g,'suppress',{'arg_direct','EEG','report_args'})));

