function [IC,TO] = Donahue_2022(acc,gyr,mag,fs)

% -------------------------------------------------------------------------
% Donahue_2022
%
% IC and TO detection from antero-posterior acceleration
%
% INPUTS:
%   acc         Nx3 acceleration [AP,ML,VT]
%   gyr         Nx3 angular velocity [AP,ML,VT]
%   mag         Nx3 magnetometer [AP,ML,VT]
%   fs          sampling frequency
%
% OUTPUT:
%   IC          detected initial contacts
%   TO          detected toe offs
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
acc_vt = acc_k(:,3);

% Compute resultant acceleration
acc_R = vecnorm(acc_k');

% Rules
sr_ic = 50;                     % spatial rule for IC detection
tr_ic = round(0.5*fs);          % temporal rule for IC detection
tr_to_lowerb = round(0.1*fs);   % lower bound for TO detection
sr_to = 29.4;                   % spatial rule for TO detection


% -------------------------------------------------------------------------

% IC detection
[~,locs] = findpeaks(acc_R,'MinPeakHeight',sr_ic,'MinPeakDistance',tr_ic);

if isempty(locs)
    IC = [IC; NaN];
    TO = [TO; NaN];
else

    % TO detection
    IC_k = locs;
        
    for k = 1:length(locs)

        lowerb = IC_k(k)+tr_to_lowerb;

        if k < length(locs)
            stride = IC_k(k+1)-IC_k(k);
        else
            stride = IC_k(k)-IC_k(k-1);  % approximate last stride
        end

        half_stride = round(stride/2);
        upperb = lowerb+half_stride;

        LB = max(lowerb,1);
        UB = min(upperb,length(acc_R));

        if LB >= UB
            TO = [TO; NaN];
        else
            segm_vt = acc_vt(LB:UB);

            % spatial rule 1
            [pks_to,locs_to] = findpeaks(segm_vt);

            if ~isempty(locs_to)
                [~,idx] = max(pks_to);
                TO_k = locs_to(idx)+(LB-1);
            else
                % spatial rule 2
                idx_thr = find(segm_vt > sr_to,1,'first');
                if ~isempty(idx_thr)
                    TO_k = idx_thr+(LB-1);
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

