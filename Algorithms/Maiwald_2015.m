function [IC] = Maiwald_2015(acc,gyr,fs)

% -------------------------------------------------------------------------
% Maiwald_2015
% IC detection from vertical acceleration using peak detection
% Windowed in mid-swing to mid-swing cycles
%
% INPUTS:
%   acc         Nx3 acceleration [VT,ML,AP]
%   gyr         Nx3 angular velocity [VT,ML,AP]
%   fs          sampling frequency
%   MS          midswing indices (vector)
%
% OUTPUT:
%   IC          detected initial contacts
% -------------------------------------------------------------------------

IC = [];

% Extract acceleration components
acc_vt = acc(:,1);
gyr_ml = gyr(:,3);

% Filtering
[b_l,a_l] = butter(2,2/(fs/2),'low');
gyr_filt = filtfilt(b_l,a_l,gyr_ml);
[~,tn] = findpeaks(gyr_filt);

[b_h,a_h] = butter(4,80/(fs/2),'high');
acc_filt = filtfilt(b_h,a_h,acc_vt);


% -------------------------------------------------------------------------

% Loop over tn to find IC
for i = 1:length(tn)

    if i == 1
        start_idx = 1;
    else
        start_idx = tn(i-1);
    end
    end_idx = tn(i);

    segment = acc_filt(start_idx:end_idx);
    [pks,locs] = findpeaks(segment);

    if isempty(pks)
        IC = [IC; NaN];
    else
        [~,idx] = max(pks);
        IC = [IC; locs(idx)+start_idx-1];
    end
    
end

IC = IC(:);

end

