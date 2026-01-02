function [IC,TO] = Mo_2018(acc_foot,acc_tibia,MS)

% -------------------------------------------------------------------------
% This implementation is based on the methodology described in Mo, 2018.
% The original authors of the referenced papers are not responsible for this implementation.
%
% IC and TO detection from acceleration using two sensors 
% (dorsum of the foot + medial tibia)
% Windowed in mid-swing to mid-swing cycles
%
% INPUTS:
%   acc_foot    Nx3 acceleration [AP,ML,VT]
%   acc_tibia   Nx3 acceleration [VT,AP,VML]
%   fs          sampling frequency
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

% Gravity normalization
g = 9.81;
acc_foot = acc_foot./g;
acc_vt = acc_tibia(:,1)./g;

% Compute resultant acceleration
acc_R = vecnorm(acc_foot');  


% -------------------------------------------------------------------------

% Create window start and end indices
win_start = [1;MS(:)];
win_end = [MS(:);length(acc_R)];


% Loop over windows
for w = 1:length(win_start)

    % IC detection
    segm_R = acc_R(win_start(w):win_end(w));
    [pks,locs] = findpeaks(segm_R);

    if isempty(locs)
        IC = [IC; NaN];
        TO = [TO; NaN];
    else
        [~,idx] = max(pks);
        IC_k = locs(idx)+win_start(w)-1;
        IC = [IC; IC_k];
    
        % TO detection
        segm_vt = acc_vt(IC_k:win_end(w));
        [~,locs_to] = findpeaks(segm_vt,'NPeaks',3);

        if length(locs_to) < 3
            TO_k = NaN;
        else
            sign_int = segm_vt(locs_to(2):locs_to(3));
            [~, locs_min] = min(sign_int);
    
            if ~isempty(locs_min)
                TO_k = locs_min+(IC_k-1)+(locs_to(2)-1);
            else
                TO_k = NaN;
            end
        end 
        TO = [TO; TO_k];
    end
end

IC = IC(:);
TO = TO(:);

end


