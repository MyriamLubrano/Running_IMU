function [IC,TO] = Reenalda_2019(acc,gyr,MS)

% -------------------------------------------------------------------------
% Reenalda_2019
%
% IC and TO detection from resultant acceleration and angular velocity
% Windowed in mid-swing to mid-swing cycles
%
% INPUTS:
%   acc         Nx3 acceleration [AP,ML,VT]
%   gyr         Nx3 angular velocity [AP,ML,VT]
%   MS          midswing indices (vector)
%
% OUTPUT:
%   IC          detected initial contacts
%   TO          detected toe offs
% -------------------------------------------------------------------------

IC = [];
TO = [];

% Extract acceleration components
acc_ap = -acc(:,1);
acc_ml = acc(:,2);
acc_vt = acc(:,3);

gyr_ap = -gyr(:,1);
gyr_ml = gyr(:,2);
gyr_vt = gyr(:,3);

% Compute resultant acceleration
acc_R = sqrt(acc_ap.^2+acc_ml.^2+acc_vt.^2);
gyr_R = sqrt(gyr_ap.^2+gyr_ml.^2+gyr_vt.^2);


% -------------------------------------------------------------------------

% Create window start and end indices
win_start = [1;MS(:)];
win_end = [MS(:);length(acc_R)];


% Loop over windows
for w = 1:length(win_start)
    
    segment_acc = acc_R(win_start(w):win_end(w));
    [pks_acc,locs_acc] = findpeaks(segment_acc);

    if isempty(pks_acc)
        IC = [IC; NaN];
        TO = [TO; NaN];
    else
        % Peak on acc magnitude: IC
        [~,idx_acc] = max(pks_acc);
        IC_k = locs_acc(idx_acc)+win_start(w)-1;
        IC = [IC; IC_k];
        
        % Peak on gyr magnitude: TO
        start_TO = locs_acc(idx_acc);
        segment_gyr = gyr_R(win_start(w)+start_TO:win_end(w));
        [pks_gyr,locs_gyr] = findpeaks(segment_gyr);
        
        if isempty(pks_gyr)
            TO = [TO; NaN];
        else
            [~,idx_gyr] = max(pks_gyr);
            TO_k = locs_gyr(idx_gyr)+win_start(w)+start_TO-1;
            TO = [TO; TO_k];
        end
    end
end

IC = IC(:);
TO = TO(:);

end

