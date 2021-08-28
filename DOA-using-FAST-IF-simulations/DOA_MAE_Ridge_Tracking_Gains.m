% close all;
clear all;
addpath('C:\Program Files\MATLAB\R2021a\toolbox\tfsa_5_5');

%% Parameters:
N_sources = 2;
N_sensors=2;
SNR=50;
Gain_range_dB = 0:-1:-20;
Gain_range = 10.^(Gain_range_dB/20);
% samples of sig
N_samples = 128;
n=0:N_samples-1;
% set Degrees parameters
theta = [-15,15]*pi/180;
theta_deg=round(theta *180/pi);
theta_range = -90:1:90;
theta_axis=zeros(1,length(theta_range));
theta_axis(theta_deg+ceil(length(theta_range)/2))=1;
% Parameters for signals:
init_freqs = [0.3 , 0.25];
curve_slope = [0.2 , -0.1];
sel_samples = 1:N_samples;

%% Create signal sources
delta = 3;
N_scenarios=10;
for gain_idx = 1:length(Gain_range_dB) 
    s(1,:) = exp(2*pi*1i*(init_freqs(1)*n+curve_slope(1)*n.^3/(N_samples*N_samples*3))).';
    s(2,:) = Gain_range(gain_idx)*exp(2*pi*1i*(init_freqs(2)*n+curve_slope(2)*n.^3/(N_samples*N_samples*3))).';
    IF_O(1,:) = init_freqs(1)+curve_slope(1)*3*n.^2/(N_samples*N_samples*3);
    IF_O(2,:) = init_freqs(2)+curve_slope(2)*3*n.^2/(N_samples*N_samples*3);

    % Channel matrix A 
    A = exp(1j*pi*[0:N_sensors-1].'*sin(theta));                
    X = A*s;                             % mixed source
    for curr_scenario=1:N_scenarios 
        % AWGN     
        sigma = 10^(-SNR/20);
        w = sigma*(randn(N_sensors,length(n)) + 1j*(randn(N_sensors,length(n))))/sqrt(2); % noise       
        X=X+w;     

        [Ridges_IF, Ridge_tracking_signals]= multi_sensor_source_separation_ridge_tracking_m(X, N_sources, delta,N_sensors);

        for curr_src=1:N_sources
            for curr_sensor=1:N_sensors
                Signals_to_DOA(curr_sensor,:)=Ridge_tracking_signals(curr_sensor,curr_src,:);
            end
            
            p=TMMUSIC(cov(Signals_to_DOA'), 2, N_sensors, 1, 1, theta_range');
            [x,y]=max(p);
            y1(curr_src)=y;
        end
        y1=y1-90;
        abs_err = abs((sort(y1)-sort(theta_deg)));
        MEAN_ABS_ERR_FAST_IF(1,curr_scenario)=abs_err(1);
        MEAN_ABS_ERR_FAST_IF(2,curr_scenario)=abs_err(2);
        
    end
    MAE_FAST_IF_vec(1,gain_idx)=mean(MEAN_ABS_ERR_FAST_IF(1,:))
    MAE_FAST_IF_vec(2,gain_idx)=mean(MEAN_ABS_ERR_FAST_IF(2,:))
end

figure;
plot(-Gain_range_dB,MAE_FAST_IF_vec(1,:),'linewidth',2); grid; hold on;
plot(-Gain_range_dB,MAE_FAST_IF_vec(2,:),'linewidth',2); hold off;
title({['Ridge Tracking Mean Absolute Error Vs difference in gain'];['SNR = ',num2str(SNR),'dB']}); legend('1st Signal','2nd Signal'); 
xlabel('Gain difference [dB]'); ylabel('Mean Absolute Error [degrees]');

