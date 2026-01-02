function [IC] = Aubol_2020(acc,fs)

% -------------------------------------------------------------------------
% This implementation is based on the methodology described in Aubol, 2020.
% The original authors of the referenced papers are not responsible for this implementation.
% 
% IC detection from resultant acceleration using peak detection
%
% INPUTS:
%   acc         Nx3 acceleration [VT,AP,ML]
%   fs          sampling frequency
%
% OUTPUT:
%   IC          detected initial contacts
%
% Author: Myriam Lubrano
% -------------------------------------------------------------------------

IC = [];

% Filtering
[b,a] = butter(4,70/(fs/2),'low');
acc_filt = filtfilt(b,a,acc);

% Compute resultant acceleration
acc_R = vecnorm(acc_filt');  

% Detect peaks in the resultant acceleration
[pk_acc,loc_acc] = findpeaks(acc_R,'MinPeakHeight',60,'MinPeakDistance',300);

% Compute resultant jerk
jerk_res = diff(acc_res);


% -------------------------------------------------------------------------

window = round(0.150*fs);       % 150 ms window before each acceleration peak

for p = 1:length(loc_acc)
    if loc_acc(p) > window
        sig_j = jerk_res(loc_acc(p)-window:loc_acc(p));
        [~, max_idx] = max(sig_j);
        ind_max_j(p) = max_idx+loc_acc(p)-window;
    else
        sig_j = jerk_res(1:loc_acc(p));
        [~, max_idx] = max(sig_j);
        ind_max_j(p) = max_idx+loc_acc(p);
    end
    
end 

window = round(0.075*fs);       % 75 ms window before each jerk peak

for p = 1:length(loc_acc)
    if ind_max_j(p) > window
        sig_a = acc_res(ind_max_j(p)-window:ind_max_j(p));
        [~, min_idx] = min(sig_a);
        ind_min_a(p) = min_idx+ind_max_j(p)-window;
    else
        sig_a = acc_res(1:ind_max_j(p));
        [~, min_idx] = min(sig_a);
        ind_min_a(p) = min_idx+ind_max_j(p);
    end
end

% Determine prominence threshold
sorted_peaks = sort(pk_acc,'descend');
if length(sorted_peaks) >= 3
    prom_th = 0.2*sorted_peaks(3);
else
    prom_th = 0.2*sorted_peaks(end);
end

contact_locs = [];
prev_jerk_loc = NaN;

% Filter local minima based on frequency and prominence criteria
for p = 1:length(ind_min_a)
    min_loc = ind_min_a(p);
    jerk_peak_loc = ind_max_j(p);
    
    if ~isnan(prev_jerk_loc)
        dt = (jerk_peak_loc-prev_jerk_loc)/fs;
        f_peak = 1/dt;
        if f_peak > 40
            continue
        end
    end
    prev_jerk_loc = jerk_peak_loc;

    if p <= length(pk_acc)
        prom = pk_acc(p)-acc_res(min_loc);
        if prom < prom_th
            continue  
        end
    end
    contact_locs(end+1) = min_loc;
end

% Select earliest local minimum prior to each acceleration peak
IC = zeros(size(loc_acc));
for i = 1:length(loc_acc)
    prior_mins = contact_locs(contact_locs < loc_acc(i));
    if ~isempty(prior_mins)
        IC(i) = prior_mins(end);  % closest minimum before peak
    else
        IC(i) = NaN;
    end
end

IC = IC(:);

end



