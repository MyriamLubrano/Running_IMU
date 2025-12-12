function [IC,TO] = Lee_2010(acc,fs)

% -------------------------------------------------------------------------
% Lee_2010
%
% IC and TO detection from antero-posterior acceleration
%
% INPUTS:
%   acc         Nx3 acceleration [VT,ML,AP]
%   fs          sampling frequency
%
% OUTPUT:
%   IC          detected initial contacts
%   TO          detected toe offs
% -------------------------------------------------------------------------

IC = [];
TO = [];

% Extract acceleration components
acc_ap = acc(:,3);

dist1 = round(0.08*fs);       % distance in samples from IC, to be adapted
dist2 = round(0.05*fs);       % distance in samples from next IC, to be adapted


% -------------------------------------------------------------------------

% IC detection
[~,locs] = findpeaks(acc_ap,'MinPeakHeight',15,'MinPeakDistance',150);     % to be adapted

if isempty(locs)
    IC = [IC; NaN];
    TO = [TO; NaN];
else
    IC = locs;
    
    win_start = locs(:)+dist1;
    win_end = [locs(2:end)-dist2;length(acc_ap)];
    
    for w = 1:length(win_start)

        start_idx = max(win_start(w),1);
        end_idx   = min(win_end(w),length(acc_ap));

        if start_idx >= end_idx
            TO = [TO; NaN];
            continue
        end
        
        window_to = acc_ap(start_idx:end_idx);
        [pks_to,locs_to] = findpeaks(window_to);
    
        % TO detection
        if isempty(locs_to)
            TO = [TO; NaN];
        else
            [~,idx_to] = max(pks_to);
            TO = [TO; locs_to(idx_to)+(win_start(w)-1)];
        end
    end
end

IC = IC(:);
TO = TO(:);

end

