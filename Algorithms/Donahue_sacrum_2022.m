function [IC,TO] = Donahue_sacrum_2022(acc,gyr,mag,fs)

% -------------------------------------------------------------------------
% This implementation is based on the methodology described in Donahue, 2022.
% The original authors of the referenced papers are not responsible for this implementation.
%
% IC and TO detection from antero-posterior acceleration
%
% INPUTS:
%   acc         Nx3 acceleration [VT,ML,AP]
%   gyr         Nx3 angular velocity [VT,ML,AP]
%   mag         Nx3 magnetometer [VT,ML,AP]
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

% Filtering
[b,a] = butter(4,35/(fs/2),'low');
acc = filtfilt(b,a,acc); 
gyr = filtfilt(b,a,gyr);
mag = filtfilt(b,a,mag);

% Kalman filter
fuse = ahrsfilter('SampleRate',fs,'ReferenceFrame','NED','OrientationFormat','Rotation matrix');
[orientation,~] = fuse(acc,gyr,mag);
R = orientation;
for i = 1:length(gyr)
    % from local to global reference system
    acc_k(i,:) = R(:,:,i)'*acc(i,:)';
    gyr_k(i,:) = R(:,:,i)'*gyr(i,:)';
end

% Extract acceleration components
acc_ap = acc_k(:,3);

% Rules
sr_ic = 5;                      % spatial rule for IC detection
tr_ic = round(0.2*fs);          % temporal rule for IC detection
tr_to_lowerb = round(0.1*fs);   % lower bound for TO detection
tr_to_upperb = round(0.02*fs);  % lower bound for TO detection


% -------------------------------------------------------------------------

% IC detection
[~,locs] = findpeaks(acc_ap,'MinPeakHeight',sr_ic,'MinPeakDistance',tr_ic);

if isempty(locs)
    IC = [IC; NaN];
    TO = [TO; NaN];
else

    % TO detection
    IC_k = locs(:);
    IC = [IC; IC_k];
    prev_range = [];
        
    for k = 1:length(locs)

        LB = IC_k(k)+tr_to_lowerb;

        if k < length(locs)
            UB = IC_k(k+1)-tr_to_upperb;
            if UB > LB
                prev_range = UB-LB;    % store window length
            end
        else                           % last stride  
            upperb = LB+prev_range;      
            UB = min(upperb,length(acc_ap));    
        end

        if LB >= UB
            TO = [TO; NaN];
        else
            segm_ap = acc_ap(LB:UB);

            % spatial rule 1
            [pks_to,locs_to] = findpeaks(segm_ap);

            if ~isempty(locs_to)
                [~,idx] = max(pks_to);
                TO_k = locs_to(idx)+(LB-1);
            else
                % spatial rule 2
                slope = diff(segm_ap);
                [~,idx_slope] = max(slope);
                
                if ~isempty(idx_slope)
                    TO_k = idx_slope+(LB-1);
                else
                    TO_k = NaN;
                end
            end
            TO = [TO; TO_k];
        end
    end
end

IC = IC(:);
TO = TO(:);

end


