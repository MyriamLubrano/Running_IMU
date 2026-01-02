function [IC,TO] = Yang_2022(acc,gyr,fs,MS)

% -------------------------------------------------------------------------
% This implementation is based on the methodology described in Yang, 2022.
% The original authors of the referenced papers are not responsible for this implementation.
%
% IC and TO detection using peak detection on both acceleration and angular
% velocity
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
%
% Authors: Myriam Lubrano, Rachele Rossanigo, Elena Dipalma 
% -------------------------------------------------------------------------

IC = [];
TO = [];

% Extract acceleration components
acc_vt = acc(:,1);
acc_ap = -acc(:,2);
acc_ml = acc(:,3);
gyr_ml = gyr(:,3);

% Compute resultant acceleration
acc_R = sqrt(acc_ap.^2+acc_ml.^2+acc_vt.^2);


% -------------------------------------------------------------------------

% Create window start and end indices
win_start = [1;MS(:)];
win_end = [MS(:);length(acc_R)];


% Loop over windows
for w = 1:length(win_start)

    int = win_start(w):win_end(w);
    segm_acc = acc_R(int);
    segm_acc = segm_acc(1:floor(length(segm_acc)/2));    
    
    if length(int) < 0.3*fs  % skip too short windows
        continue
    end

    % IC detection
    SDA = [];                       % square difference of acceleration
    for i = 1:length(segm_acc)-1
        SDA(end+1) = segm_acc(i+1)-segm_acc(i);
    end
    SDA = SDA.^2;
    
    [pks,locs] = findpeaks(SDA);   

    if isempty(locs)
        [~,IC_k] = max(SDA);        % fallback
    else
        [~,idx] = max(pks);
        IC_k = locs(idx);   
    end
    IC = [IC; IC_k+win_start(w)-1];

    % TO detection
    int_to = int(IC_k):int(end);
    segm_gyr = -gyr_ml(int_to);

    % midswing 
    cutoff = 30;
    Wn = cutoff/(fs/2);
    [B,A] = butter(2,Wn,'low'); 

    gyr_filt = filtfilt(B,A,segm_gyr);
    [~, locs_to] = findpeaks(gyr_filt);

    if length(locs_to) == 1
        TO_k = locs_to;
    elseif length(locs_to) > 1
        TO_k = locs_to(2);
    elseif isempty(locs_to)
        TO_k = NaN;
    end
    TO_k = TO_k+IC_k;

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




