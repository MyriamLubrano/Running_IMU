function [IC,TO] = Chew_2018(acc,MS)

% -------------------------------------------------------------------------
% This implementation is based on the methodology described in Chew, 2018.
% The original authors of the referenced papers are not responsible for this implementation.
%
% IC and TO detection from antero-posterior acceleration 
% Windowed in mid-swing to mid-swing cycles
%
% INPUTS:
%   acc         Nx3 acceleration from foot sensor [AP,ML,VT]
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

% Extract antero-posterior acceleration
acc_ap = -acc(:,1);


% -------------------------------------------------------------------------

% Create window start and end indices
win_start = [1;MS(:)];
win_end = [MS(:);length(acc_ap)];


% Loop over windows
for w = 1:length(win_start)
    
    segment = -acc_ap(win_start(w):win_end(w));
    [pks,locs] = findpeaks(segment);

    if isempty(pks)
        IC = [IC; NaN];
        TO = [TO; NaN];
    elseif length(pks) == 1
        IC = [IC; locs(1)+win_start(w)-1];
        TO = [TO; NaN];
    else
        % Take first larger peak as IC
        [~,idx] = max(pks);
        IC = [IC; locs(idx)+win_start(w)-1];
        
        % Second peak as TO
        [~, idx2] = max(pks(idx+1:end));
        TO = [TO; locs(idx+idx2)+win_start(w)-1];
    end
end

IC = IC(:);
TO = TO(:);

end


