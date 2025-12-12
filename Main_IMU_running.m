%% Main: Temporal Parameters Estimation
% Example of main script for temporal parameters estimation using IMUs

clc, clear, close all


%% Loading data

fs = 800;       % sampling frequency [Hz]

acc = ;         % acceleration data (Nx3) [m/s2] 
gyr = ;         % angular velocity data (Nx3) [rad/s]


%% Segmentation into midswing to midswing cycles

gyr_ml = gyr(:,2);      % medio-lateral angular velocity
[r,lags] = xcorr(gyr_ml,'normalized');
[pks,locs] = findpeaks(r,'MinPeakHeight',0); 
t = mean(diff(locs/fs));
strFr = 1/t; 
cutoff = 0.6*strFr;     % adaptive filter

Wn = cutoff/(fs/2);
[b,a] = butter(2,Wn,'low'); 
gyr_filt = filtfilt(b,a,gyr_ml);
[~,MS] = findpeaks(gyr_filt,'MinPeakHeight',0,'MinPeakDistance',30);


%% Algorithm implementation: IC and TO detection
% acc and gyr should go from 1 to length(acc), with MS representing the
% indices of midswing instants inside the interval

% 'algorithm' is a placeholder for the chosen function
% e.g., Blauberger_2021, Chew_2018, etc.
[IC_samples,TO_samples] = algorithm(acc,gyr,fs,MS);

% convert sample indices to seconds
IC = IC_samples/fs;     
TO = TO_samples/fs;


%% Temporal parameters estimation
% Check the length of the vectors

% stride duration: difference between consecutive ICs [s]
stride = (IC(2:end)-IC(1:end-1));   % one side
stride = (IC(3:end)-IC(1:end-2));   % both sides
% stance duration: difference between TO and IC of the same step [s]
stance = TO-IC;

