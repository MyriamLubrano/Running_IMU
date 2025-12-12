function [IC,TO] = Falbriard_2018(acc,gyr,fs)

% -------------------------------------------------------------------------
% Falbriard_2018
%
% IC and TO detection from resultant acceleration
% Windowed in mid-swing to mid-swing cycles Falbriard 2018
%
% INPUTS:
%   acc         Nx3 acceleration [AP,ML,VT]
%   fs          sampling frequency
%
% OUTPUT:
%   IC          detected initial contacts
%   TO          detected toe offs
% -------------------------------------------------------------------------

IC = [];
TO = [];

% Midswing segmentation
gyr_ml = gyr(:,2);
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
acc = filtfilt(b,a,acc);

% Extract acceleration components
acc_ap = -acc(:,1);
acc_ml = acc(:,2);
acc_vt = acc(:,3);

% Compute resultant acceleration
acc_R = sqrt(acc_ap.^2+acc_ml.^2+acc_vt.^2);


% -------------------------------------------------------------------------

% Create window start and end indices
win_start = MS(1:end-1);    
win_end   = MS(2:end);    


% Loop over windows
for w = 1:length(win_start)
    
    segment = acc_R(win_start(w):win_end(w));
    [pks,locs] = findpeaks(segment);

    if isempty(pks)
        IC = [IC; NaN];
        TO = [TO; NaN];
    elseif length(pks) == 1
        IC = [IC; locs(1)+win_start(w)-1];
        TO = [TO; NaN];
    else
        % Take first peak as IC, second peak as TO
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

