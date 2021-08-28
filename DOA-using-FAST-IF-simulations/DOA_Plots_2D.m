% close all;
clear all; 
% close all;
% Need to have to following tool box:
addpath('C:\Program Files\MATLAB\R2021a\toolbox\tfsa_5_5');

%% Parameters:
N_sensors=4;
N_sources = 4;
rand_sources = 1;
% samples of sig
Ts = 1;
fs = 1/Ts;
time_axis = 0:Ts:128-Ts;
N_samples = length(time_axis);
n=0:N_samples-1;
nT = n*Ts;
% generate noise
SNR=0;
sigma = 10^(-SNR/20);
% Degrees parameters
angle_range = 0:1:180;
% define angles
theta_vec = [30, 140];
phi_vec = [45 , 140];

%% Create signal sources
% [s, IF_O, N_sources] = Create_FM_Signals(3, n, Ts, N_sources);    % 1 - original, 2 - cyclic chirp

%% calculate channel and antennas signals
N_Scenarios = 100;
plot_flag = 0;
for curr_scenario = 1:N_Scenarios
    if rand_sources
        new_N_sources = randi(N_sources);
    else
        new_N_sources = N_sources;
    end
    % Create new sources
    [s, IF_O, N_sources] = Create_FM_Signals(3, n, Ts, new_N_sources);    % 1 - original, 2 - cyclic chirp
    phi_rad = phi_vec(randi(length(phi_vec),[1,N_sources]))*pi/180;
    phi_deg = round(phi_rad *180/pi);
    theta_rad = theta_vec(randi(length(theta_vec),[1,N_sources]))*pi/180;
    theta_deg = round(theta_rad *180/pi);
    phi_axis=zeros(1,length(angle_range));
    theta_axis=zeros(1,length(angle_range));
    phi_axis(phi_deg+1) = 1;
    theta_axis(theta_deg+1) = 1;

    % Channel matrix for plain Z
    A_z = exp(1j*pi*[0:N_sensors-1].'*cos(theta_rad));  
    Z_Sig = A_z*s.';

    % Additive WGN
    w = sigma*(randn(N_sensors,length(n)) + 1j*(randn(N_sensors,length(n))))/sqrt(2); % noise
    Z_Sig=Z_Sig+w;

    % Channel matrix for plain XY
    A_x = exp(1j*pi*[0:N_sensors-1].'*(cos(phi_rad).*sin(theta_rad)));  
    X_Sig = A_x*s.';

    % Additive WGN
    w = sigma*(randn(N_sensors,length(n)) + 1j*(randn(N_sensors,length(n))))/sqrt(2); % noise
    X_Sig=X_Sig+w;

    %% Fast-IF and DOA estimation
    win_length = 33;
    delta = 6;
    L = 100;
    relative_energy_th = 0.1; % zero if we don't want to estimate number of sources
    [theta_est,phi_est, N_Sources_est] = FAST_IF_2D_DOA(Z_Sig,X_Sig, N_sensors,N_sources,win_length,delta,L,relative_energy_th,plot_flag,[theta_deg,phi_deg]);
    
    % figure of DOA
    figure(100);
    colors_vec = {'blue', 'green', 'black'};
    for curr_source = 1:N_sources
        scatter(theta_est(curr_source),phi_est(curr_source),50,'red','x'); hold on;
        scatter(theta_deg(curr_source),phi_deg(curr_source),50,colors_vec{curr_source},'o'); hold on;
    end
    title({'2D FAST-IF DOA';...
           ['SNR = ',num2str(SNR),'dB']}); xlabel('\theta [Degrees]'); ylabel('\phi [Degrees]'); grid on;
    xlim([0 180]);
    ylim([0 180]);
end



