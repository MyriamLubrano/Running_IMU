function [IC] = VanWerkhoven_2019(acc,MS)

% -------------------------------------------------------------------------
% VanWerkhoven_2019
% IC detection from resultant acceleration using peak detection
% Windowed mid-swing → mid-swing search
%
% INPUTS:
%   acc     Nx3 acceleration [AP,ML,VT]
%   MS          midswing indices (vector)
%
% OUTPUT:
%   IC          detected initial contacts
% -------------------------------------------------------------------------

IC = [];

% Extract acceleration components
acc_ap = -acc(:,1);
acc_ml = acc(:,2);
acc_vt = acc(:,3);

% Compute resultant acceleration
acc_R = sqrt(acc_ap.^2+acc_ml.^2+acc_vt.^2);


% -------------------------------------------------------------------------

% Create window start and end indices
win_start = [1;MS(:)];
win_end = [MS(:);length(acc_R)];


% Loop over windows
for w = 1:length(win_start)
    
    segment = acc_R(win_start(w):win_end(w));
    [pks,locs] = findpeaks(segment);

    if isempty(pks)
        IC = [IC; NaN];
    else
        [~,idx] = max(pks);
        IC = [IC; locs(idx)+win_start(w)-1];
    end
end

IC = IC(:);

end

