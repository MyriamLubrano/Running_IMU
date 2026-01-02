function [IC,TO] = Falbriard_2020(gyr,fs)

% -------------------------------------------------------------------------
% This implementation is based on the methodology described in Falbriard, 2020.
% The original authors of the referenced papers are not responsible for this implementation.
%
% IC and TO detection from resultant angular velocity
% Windowed in mid-swing to mid-swing cycles as in Falbriard 2018
%
% INPUTS:
%   gyr         Nx3 angular velocity [AP,ML,VT]
%   fs          sampling frequency
%
% OUTPUT:
%   IC          detected initial contacts
%   TO          detected toe offs
%
% Author: Myriam Lubrano
% -------------------------------------------------------------------------

IC = [];
TO = [];

% Extract acceleration components
gyr_ml = gyr(:,2);

% Midswing segmentation through an adaptive low-pass filter
stride_autocorr = xcorr(gyr_ml,'coeff');
stride_autocorr = stride_autocorr(length(gyr_ml):end);   % positive lags

[~,peaks] = findpeaks(stride_autocorr,'MinPeakDistance',round(fs*0.5)/2);       % one side
stride_frequency = fs/mean(diff(peaks)); 

cutoff_freq = 0.6*stride_frequency;

[b,a] = butter(2,cutoff_freq/(fs/2),'low');
gyr_ml_filt = filtfilt(b,a,gyr_ml);
[~,loc_MS] = findpeaks(gyr_ml_filt);     % loc_MS: midswing indices [samples]
MS = loc_MS;                        

% Filtering
[b,a] = butter(2,30/(fs/2),"low");
gyr_ml = filtfilt(b,a,gyr_ml);


% -------------------------------------------------------------------------

% Create window start and end indices
win_start = MS(1:end-1);
win_end   = MS(2:end);


% Loop over windows
for w = 1:length(win_start)
    
    segment = -gyr_ml(win_start(w):win_end(w));
    [pks,locs] = findpeaks(segment);

    if isempty(locs) || length(locs) < 2
        IC = [IC; NaN];
        TO = [TO; NaN];
    else
        % Take first peak as IC
        [~,idx] = max(pks);
        IC = [IC; locs(idx)+win_start(w)-1];
        
        % Second peak as TO
        [~, idx2] = max(pks(idx+1:end));
        TO = [TO; locs(idx+idx2)+win_start(w)-1];
    end
end

IC = IC(:);
TO = TO(:);

end


