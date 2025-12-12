function [IC,TO] = Mitschke_2018(acc,gyr,fs,MS)

% -------------------------------------------------------------------------
% Mitschke_2018
%
% IC and TO detection from acceleration and angular velocity
% Windowed in mid-swing to mid-swing cycles
%
% INPUTS:
%   acc         Nx3 acceleration [VT,AP,ML]
%   gyr         Nx3 angular velocity [VT,AP,ML]
%   fs          sampling frequency
%   MS          midswing indices (vector)
%
% OUTPUT:
%   IC          detected initial contacts
%   TO          detected toe offs
% -------------------------------------------------------------------------

IC = [];
TO = [];

% Filtering
[b,a] = butter(4,200/(fs/2),'low');
acc = filtfilt(b,a,acc);

[b,a] = butter(4,50/(fs/2),'low');
gyr = filtfilt(b,a,gyr);

% Gravity normalization
g = 9.81;
acc = acc./g;

% Extract acceleration components
acc_ap = -acc(:,2);
acc_vt = acc(:,1);
gyr_ap = -gyr(:,2);

% High-pass filter for vertical acceleration
[b,a] = butter(4,80/(fs/2),'high');
acc_vt = filtfilt(b,a,acc_vt);


% -------------------------------------------------------------------------

% Create window start and end indices
win_start = [1;MS(:)];
win_end = [MS(:);length(acc_vt)];


% Loop over windows
for w = 1:length(win_start)

    % IC detection
    segm_acc_vt = acc_vt(win_start(w):win_end(w));
    [pks,locs] = findpeaks(segm_acc_vt);

    if isempty(locs)
        IC = [IC; NaN];
    else
        [~,idx] = max(pks);
        IC = [IC; locs(idx)+win_start(w)-1];
    end

    % TO detection
    segm_acc_ap = acc_ap(win_start(w):win_end(w));
    segm_gyr_ap = gyr_ap(win_start(w):win_end(w));

    acc_hor = zeros(1,length(segm_acc_ap));

    int = locs(idx):length(segm_acc_ap);
    int_gyr = segm_gyr_ap(int); 
    int_gyr = int_gyr-mean(int_gyr);
    t_int = int/fs;
    
    theta = cumtrapz(t_int,int_gyr)';
    theta = theta-(theta(end)-theta(1));
    
    acc_hor(int) = segm_acc_ap(int).*cos(theta')-segm_acc_vt(int).*sin(theta');

    % Apply 2g threshold
    th = 2; 
    dist = fs*0.06;      % to avoid peaks too close to the main one and to be adapted
    int = locs(idx)+dist:length(acc_hor);
    ind_th = acc_hor(int) >= th;
    pass = find(diff(ind_th) == 1)';
    
    if isempty(pass) 
        TO = [TO; NaN];
    else
        TO = [TO; pass(1)+(locs(idx)+dist-1)+(win_start(w)-1)];  
    end
end

IC = IC(:);
TO = TO(:);

end

