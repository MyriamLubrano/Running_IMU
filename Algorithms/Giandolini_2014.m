function [IC_min] = Giandolini_2014(acc_heel,acc_met,fs,MS)

% -------------------------------------------------------------------------
% Giandolini_2014
%
% IC detection from vertical acceleration of heel and forefoot
% Windowed in mid-swing to mid-swing cycles
%
% INPUTS:
%   acc_heel    Nx3 acceleration from heel sensor [VT,AP,ML]
%   acc_met     Nx3 acceleration from metatarsal head sensor [VT,AP,ML]
%   fs          sampling frequency
%   MS          midswing indices (vector)
%
% OUTPUT:
%   IC          detected initial contacts
% -------------------------------------------------------------------------

IC = [];

% Extract vertical acceleration
acc_vt_H = acc_heel(:,1);
acc_vt_M = acc_met(:,1);

% Filtering
[b,a] = butter(2,50/(fs/2),'low');
acc_H = filtfilt(b,a,acc_vt_H);
acc_M = filtfilt(b,a,acc_vt_M);


% -------------------------------------------------------------------------

% Create window start and end indices
win_start = [1;MS(:)];
win_end = [MS(:);length(acc_H)];

% Loop over windows
for w = 1:length(win_start)
    
    segment_H = acc_H(win_start(w):win_end(w));
    segment_M = acc_M(win_start(w):win_end(w));

    % Heel IC detection
    [pks_H, locs_H] = findpeaks(segment_H);
    if isempty(pks_H)
        IC_H = NaN;
    else
        [~, idx] = max(pks_H);
        IC_H = locs_H(idx)+win_start(w)-1;
    end

    % Metatarsal head IC detection
    [pks_M, locs_M] = findpeaks(segment_M);
    if isempty(pks_M)
        IC_M = NaN;
    else
        [~, idx] = max(pks_M);
        IC_M = locs_M(idx)+win_start(w)-1;
    end

    % Earliest IC
    IC = [IC; min(IC_H,IC_M)];
end

IC = IC(:);

end

