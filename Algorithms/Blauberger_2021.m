function [IC,TO] = Blauberger_2021(acc,gyr,fs,MS,min_dist)

% -------------------------------------------------------------------------
% This implementation is based on the methodology described in Blauberger, 2021.
% The original authors of the referenced papers are not responsible for this implementation.
%
% IC and TO detection from resultant acceleration and angular velocity
% Windowed in mid-swing to mid-swing cycles
%
% INPUTS:
%   acc         Nx3 acceleration [VT,AP,ML]
%   gyr         Nx3 angular velocity [VT,AP,ML]
%   MS          midswing indices (vector)
%
% OUTPUT:
%   IC          detected initial contacts
%   TO          detected toe offs
%
% Author: Myriam Lubrano
% -------------------------------------------------------------------------

IC = [];
TO = [];

% Filtering
[b,a] = butter(2,70/(fs/2),"low");
acc = filtfilt(b,a,acc);
gyr = filtfilt(b,a,gyr);

% Extract acceleration components
acc_ap = -acc(:,2);
acc_ml = acc(:,3);
acc_vt = acc(:,1);

gyr_ap = -gyr(:,2);
gyr_ml = gyr(:,3);
gyr_vt = gyr(:,1);

% Compute resultant acceleration
acc_R = sqrt(acc_ap.^2+acc_ml.^2+acc_vt.^2);
gyr_R = sqrt(gyr_ap.^2+gyr_ml.^2+gyr_vt.^2);

% Default value
if nargin < 5 || isempty(min_dist)
    min_dist = 0.06*fs;      
end


% -------------------------------------------------------------------------

% Create window start and end indices
win_start = [1;MS(:)];
win_end = [MS(:);length(acc_R)];


% Loop over windows
for w = 1:length(win_start)
    
    seg_acc = acc_R(win_start(w):win_end(w));
    seg_gyr = gyr_R(win_start(w):win_end(w));

    % Minimum in acc: IC
    [~,loc_min] = min(seg_acc);
    IC_k = loc_min+win_start(w)-1;
    IC = [IC; IC_k];

    % Local minimum in gyr: TO
    start_TO = loc_min+min_dist;
    int_gyr = seg_gyr(start_TO:end); 
    [~,loc_p] = findpeaks(int_gyr,'MinPeakHeight',5,'MinPeakDistance',20);
    
    if length(loc_p) >= 2
        interv = loc_p(1):loc_p(2);
        [~,min_loc] = min(int_gyr(interv));
        TO_k = min_loc+(loc_p(1)-1)+(start_TO-1)+(win_start(w)-1);
    else
        TO_k = NaN;
    end
    TO = [TO; TO_k];
end

IC = IC(:);
TO = TO(:);

end


