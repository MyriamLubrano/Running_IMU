function [IC,TO] = Patoz_2022(acc,fs,body_mass)

% -------------------------------------------------------------------------
% Patoz_2022
%
% IC and TO detection from vertical acceleration
% Windowed in mid-swing to mid-swing cycles
%
% INPUTS:
%   acc         Nx3 acceleration [VT,ML,AP]
%   body_mass   body mass of the participant [kg]
%
% OUTPUT:
%   IC          detected initial contacts
%   TO          detected toe offs
% -------------------------------------------------------------------------

IC = [];
TO = [];

% Reorder acceleration components as follows:
% x: ML pointing to the subject's right
% y: AP pointing anteriorly
% z: VT pointing upward

acc_imu(:,1) = -acc(:,2);
acc_imu(:,2) = -acc(:,3);
acc_imu(:,3) = acc(:,1);


% Filtering (using a truncated Fourier series)
n = length(acc_imu);
f_tr = 0.5;
g = 9.81;
N = round(n*f_tr/fs);

sig_ifft = [];
acc_filt = [];

for comp = 1:3
    fft_acc = fft(acc_imu(:,comp));
    sig_fft_trunc = zeros(size(fft_acc));
    sig_fft_trunc(1:N) = fft_acc(1:N);

    sig_ifft = real(ifft(sig_fft_trunc));
    acc_filt(:,comp) = sig_ifft(1:n);
end

med_acc = median(acc_filt,1);
g_lab = [0,0,g];
theta_deg = acosd(dot(med_acc,g_lab)/(norm(med_acc)*norm(g_lab)));

g_target = [0,0,1];  
g_measured = med_acc/norm(med_acc);  

rotation_axis_angle = vrrotvec(g_measured,g_target); 
R = vrrotvec2mat(rotation_axis_angle);  

acc_aligned = (R*acc_imu')';
acc_detrended = acc_aligned-mean(acc_aligned);

f_tr = 5;

N = round(n*f_tr/fs);
sig_ifft = [];
acc_filt = [];

for comp = 1:3
    fft_acc = fft(acc_detrended(:,comp));
    sig_fft_trunc = zeros(size(fft_acc));
    sig_fft_trunc(1:N) = fft_acc(1:N);

    sig_ifft = real(ifft(sig_fft_trunc));
    acc_filt(:,comp) = sig_ifft(1:n);
end

% Gravity normalization
data_raw(:,1) = (0:length(acc_filt)-1)/fs*1000;  % t
data_raw(:,2) = acc_filt(:,1);                   % x AP 
data_raw(:,3) = acc_filt(:,2);                   % y VT
data_raw(:,4) = acc_filt(:,3);                   % z ML

 
% -------------------------------------------------------------------------

tseries = data_raw(:,3).*body_mass;     % [N]
th = 20;                                % [N]
binary_th = tseries>th;

IC = find(diff([0;binary_th]) == 1);   
TO = find(diff([binary_th;0]) == -1); 

IC = IC(:);
TO = TO(:);

end

