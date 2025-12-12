function [IC,TO] = Schmidt_2016(acc,gyr,fs,MS)

% -------------------------------------------------------------------------
% Schmidt_2016
%
% IC and TO detection from vertical acceleration
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

% Gravity normalization
g = 9.81;
acc = acc./g;

% Extract vertical acceleration
acc_vt = acc(:,1);

% Parameters
th_acc = 5;                     % 5g threshold
dead_time = round(0.09*fs);     % 90 ms
search_win = round(0.15*fs);    % 150 ms window
tol = 0.1;                      % tolerance for considering slope as constant, to be adapted
dist = 0.07*fs;                 % to be adapted

% Check for the gyroscope
use_gyr = (nargin == 4) && ~isempty(gyr);


% -------------------------------------------------------------------------

% Create window start and end indices
win_start = [1;MS(:)];
win_end = [MS(:);length(acc_vt)];


% Loop over windows
for w = 1:length(win_start)

    % IC detection
    segm_acc = acc_vt(win_start(w):win_end(w));
    
    [pks,locs] = findpeaks(segm_acc);
    if isempty(locs)
        IC = [IC; NaN];
        TO = [TO; NaN];
        continue
    end
    [~,idx] = max(pks);
    pk_loc = locs(idx);

    if use_gyr

        gyr_z  = gyr(:,3);
        
        % Candidate points where acceleration exceeds threshold
        cand = find(segm_acc > th_acc);
        cand = cand(cand <= pk_loc);
    
        IC_k = NaN;
    
        if ~isempty(cand)
            
            % Compute local slope of angular velocity
            segm_gyr = gyr_z(win_start(w):win_end(w));
            slope = diff(segm_gyr);
    
            % Consider monotone regions
            slope_flag = abs(slope) < tol;
            slope_flag = [slope_flag; false];  % match length of segm_acc
            
            % Intersection of conditions
            valid_idx = intersect(cand,find(slope_flag));
            
            if ~isempty(valid_idx)
                % Take first valid index as lower bound of search
                lower_bound = valid_idx(1);
            else
                % Fallback: 30 ms before peak
                lower_bound = max(pk_loc-round(0.03*fs),1);
            end
            
            % Minimum in the range
            range = lower_bound:pk_loc;
            [~,lmin] = min(seg_acc(range));
            IC_k = range(lmin);
        end    
    else
        % fallback (no gyroscope)
        range = max(pk_loc-dist,1):pk_loc;
        [~,lmin] = min(segm_acc(range));
        IC_k = range(lmin);

    end
    
    IC = [IC; IC_k+(win_start(w)-1)];

    % TO detection
    start_TO = IC(end)+dead_time;
    end_TO = start_TO+search_win;

    [~,to_local] = min(segm_acc(start_TO:end_TO));

    TO_k = to_local+start_TO-1;
    TO = [TO; TO_k+(win_start(w)-1)];
end

IC = IC(:);
TO = TO(:);

end

