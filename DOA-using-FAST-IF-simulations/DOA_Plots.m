% close all;
clear all;
% Need to have to following tool box:
addpath('C:\Program Files\MATLAB\R2021a\toolbox\tfsa_5_5');

%% Parameters:
N_sensors=3;
% samples of sig
Ts = 1;
fs = 1/Ts;
time_axis = 0:Ts:128-Ts;
N_samples = length(time_axis);
n=0:N_samples-1;
nT = n*Ts;
% generate noise
SNR=100;
% set Degrees parameters
theta = [45,15]*pi/180;
theta_deg=round(theta *180/pi);
theta_range = -90:1:90;
theta_axis=zeros(1,length(theta_range));
theta_axis(theta_deg+ceil(length(theta_range)/2))=1;

sel_samples = 1:N_samples;

%% Create signal sources
[s, IF_O, N_sources] = Create_FM_Signals(1, n, Ts);    % 1 - original, 2 - cyclic chirp

n = n(sel_samples);
% Channel matrix A
A = exp(1j*pi*[0:N_sensors-1].'*sin(theta));  
X = A*s.';

% Additive WGN
sigma = 10^(-SNR/20);
w = sigma*(randn(N_sensors,length(n)) + 1j*(randn(N_sensors,length(n))))/sqrt(2); % noise
X=X+w;

%% Fast IF algorithm:
[Fast_IF_IF,Fast_IF_signals] = Multi_Sensor_FAST_IF(X,N_sensors,65, N_sources, 6,100,0,0);

figure; 
plot(Fast_IF_IF'); hold on; plot(IF_O','color','r','LineStyle','-.','LineWidth',0.75);
title('Fast-IF Signal Estimation'); xlabel('Time [sec]'); ylabel('Frequency [Hz]');

for curr_source=1:N_sources
    for curr_sensor=1:N_sensors
        Signals_to_DOA(curr_sensor,:)=Fast_IF_signals(curr_sensor,curr_source,:);
    end
    
    p=TMMUSIC(cov(Signals_to_DOA'), 2, N_sensors, 1, 1, theta_range');
    P(curr_source,:)=p;
end

figure;
plot(theta_range,P')
hold on; stem(theta_range,theta_axis)
title('FAST-IF DOA'); xlabel('Degrees'); ylabel('Spatial Spectrum');

%% Ridge tracking algorithm:
[Ridges_IF, Ridge_tracking_signals] = multi_sensor_source_separation_ridge_tracking_m(X, N_sources, 3,N_sensors);

figure; 
plot(Ridges_IF'); hold on; plot(IF_O','color','r','LineStyle','-.','LineWidth',0.75);
title('Ridge tracking Signal Estimation'); xlabel('Time [sec]'); ylabel('Frequency [Hz]'); 

for curr_source=1:N_sources
    for curr_sensor=1:N_sensors
        Signals_to_DOA(curr_sensor,:)=Ridge_tracking_signals(curr_sensor,curr_source,:);
    end
    
    p=TMMUSIC(cov(Signals_to_DOA'), 2, N_sensors, 1, 1, theta_range');
    P(curr_source,:)=p;
end

figure;
plot(theta_range,P')
hold on; stem(theta_range,theta_axis)
title('Ridge tracking DOA'); xlabel('degrees'); ylabel('Spatial Spectrum');

%% SADTFD+Viterbi algorithm:
[SADTFD_signals,SADTFD_IF] = multi_sensor_source_separation_spatial_TF_direction(X, N_sources, 3,N_sensors);
figure; 
plot(SADTFD_IF'); hold on; plot(IF_O','color','r','LineStyle','-.','LineWidth',0.75);
title('SADTFD Signal Estimation'); xlabel('Time [sec]'); ylabel('Frequency [Hz]'); 


for curr_source=1:N_sources
    for curr_sensor=1:N_sensors
        Signals_to_DOA(curr_sensor,:)=SADTFD_signals(curr_sensor,curr_source,:);
    end
    
    p=TMMUSIC(cov(Signals_to_DOA'), 2, N_sensors, 1, 1, theta_range');
    P(curr_source,:)=p;
end

figure;
plot(theta_range,P')
hold on; stem(theta_range,theta_axis)
title('SADTFD + Viterbi DOA'); xlabel('degrees'); ylabel('Spatial Spectrum');

%% TF MUSIC DOA estimation
D   = mtfd(X, 'ckd',1, 0.05, 0.05, length(X));
% Averaged Auto-TFD
D_avg = zeros(length(X), length(X));
for mm = 1:N_sensors, D_avg = D{mm,mm} + D_avg; end
D_avg = D_avg./N_sensors;
% Selection of high-energy (t,f) points
thr = 0.4*max(max(D_avg));
Tr = abs(D_avg) >= thr;
[F_trace, ~] = find(Tr);
n_p = length(F_trace);
D_s = zeros(N_sensors, N_sensors);
for m1 = 1:N_sensors
    for m2 = 1:N_sensors
        D_s(m1,m2) = (1/n_p).*sum(sum(D{m1,m2}.*Tr));
    end
end

%P=TMMUSIC(cov(X'), 2, N_sensors, n_sources, 1, theta1');
P = tf_music(D_s, N_sources, N_sensors, 2,1, theta_range);

figure;
plot(theta_range,P')
hold on; stem(theta_range,theta_axis)
title('TF-MUSIC'); xlabel('degrees'); ylabel('Spatial Spectrum');

