function [IC,TO] = Bailey_2015(acc,gyr,fs,MS)

% -------------------------------------------------------------------------
% Bailey_2015
%
% IC and TO detection from medio-lateral angular velocity
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

% Extract acceleration components
acc_vt = acc(:,1);
gyr_ml = gyr(:,3);


% -------------------------------------------------------------------------

% Create window start and end indices
win_start = [1;MS(:)];
win_end = [MS(:);length(gyr_ml)];


% Loop over windows
for w = 1:length(win_start)

    int = win_start(w):win_end(w);
    segm_gyr = gyr_ml(int);
    segm_acc = acc_vt(int);
    
    if length(int) < 0.3*fs  % skip too short windows
        continue
    end

    % IC detection
    der = diff(segm_gyr);
    [pks,locs] = findpeaks(-der);   % minima in the window

    if isempty(locs)
        IC = [IC; NaN];
    else
        [~,idx] = max(pks);
        glob_min = locs(idx);       % global minimum in the derivative

        % backward search for threshold crossing
        loc_max = segm_gyr > -0.5;
        IC_k = [];
        for j = glob_min+1:-1:2
            if loc_max(j) == 0 && loc_max(j-1) == 1
                IC_k = j;
                break
            end
        end
        if isempty(IC_k) || IC_k > length(int)/2
            IC = [IC; NaN];
        else
            IC = [IC; IC_k+win_start(w)-1];
        end
    end
    
    % TO detection
    % Segmentation through zero-crossing
    zc = [];
    for i = 1:length(segm_gyr)-1
        if segm_gyr(i) > 0 && segm_gyr(i+1) < 0 || segm_gyr(i) < 0 && segm_gyr(i+1) > 0
            zc(end+1) = i;
        end
    end
    if isempty(zc) || length(zc) < 2      % fallback
        w_start = fix(length(int)/2);
        w_end = length(int);
    else
        c0 = zc(1);
        c1 = zc(end);
        w_start = round(c0+(c1+c0)/2);
        w_end = round(w_start+(4/5)*(c1-w_start));
    end
    if w_end > length(int)
        w_end = length(int);
    end
    wind = w_start:w_end; 
    if isempty(wind) || length(wind) < 3 || length(wind) > length(int)
        wind = fix(length(int)/2):length(int);
    end

    % find local maximum
    [pks_to,locs_to] = findpeaks(segm_acc(wind));
    if isempty(locs_to)
        [~,idx_to] = max(segm_acc(wind));      % fallback
        TO_k = wind(idx_to)+w_start-1;
    else
        [~,idx_to] = max(pks_to);
        TO_k = locs_to(idx_to)+w_start-1;
    end

    if ~isnan(IC_k) && ~isnan(TO_k)
        if TO_k-IC_k < 10/100*fs || TO_k-IC_k > 70/100*fs   % TO too close to IC
            TO = [TO; NaN];
        else
            TO = [TO; win_start(w)+TO_k-1];            
        end
    else
        TO = [TO; NaN];
    end
end

IC = IC(:);
TO = TO(:);

end



