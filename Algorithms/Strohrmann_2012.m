function [IC,TO] = Strohrmann_2012(acc,fs,MS,th)

% -------------------------------------------------------------------------
% Strohrmann_2012
%
% IC and TO detection from acceleration
% Windowed in mid-swing to mid-swing cycles
%
% INPUTS:
%   acc         Nx3 acceleration [AP,ML,VT]
%   fs          sampling frequency
%   MS          midswing indices (vector)
%   th          threshold on acceleration (i.e., 2g in Strohrmann, 2012)
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

% Compute resultant acceleration
acc_R = vecnorm(acc');  

dist = 0.08;                  % distance [s] from the IC, to be adapted
dist = round(dist*fs);        % distance in samples


% -------------------------------------------------------------------------

% Create window start and end indices
win_start = [1;MS(:)];
win_end = [MS(:);length(acc_R)];


% Loop over windows
for w = 1:length(win_start)

    % IC detection
    segment = acc_R(win_start(w):win_end(w));
    [pks,locs] = findpeaks(segment);

    if isempty(locs)
        IC = [IC; NaN];
        TO = [TO; NaN];
    else
        [~,idx] = max(pks);
        IC = [IC; locs(idx)+win_start(w)-1];

        % TO detection
        start_TO = locs(idx)+dist;
        segm_TO = segment(start_TO:end);
        
        % Apply threshold
        ind_th = find(segm_TO >= th,1,'first');
        if isempty(ind_th)
            TO = [TO; NaN];
        else
            TO = [TO; ind_th+(start_TO-1)+(win_start(w)-1)];  
        end
    end
end

IC = IC(:);
TO = TO(:);

end

