close all;
clear all;
addpath('C:\Program Files\MATLAB\R2021a\toolbox\tfsa_5_5');

%% Parameters:
Ts = 1;
N_sensors=4;
SNR=0:2:20;
% samples of sig
N_samples = 128;
n=0:N_samples-1;
% set Degrees parameters
theta = [-5,15]*pi/180;
theta_deg=round(theta *180/pi);
theta_range = -90:1:90;
theta_axis=zeros(1,length(theta_range));
theta_axis(theta_deg+ceil(length(theta_range)/2))=1;
% Parameters for signals:
init_freqs = [0.6 , 0.4];
curve_slope = [0.45 , -0.45];
sel_samples = 1:N_samples;

%% Create signal sources
N_scenarios=1000;
index=0;
delta=2;
for SNR_idx=1:length(SNR)
    curr_SNR = SNR(SNR_idx);
    for curr_scenario=1:N_scenarios 
        curr_N_sources = randi(N_sensors);
        % set Degrees parameters
        theta = (randi(171,[1,curr_N_sources])-91)*pi/180;
        theta_deg=round(theta *180/pi);
        theta_axis=zeros(1,length(theta_range));
        theta_axis(theta_deg+ceil(length(theta_range)/2))=1;
        [s, IF_O, N_sources] = Create_FM_Signals(3, n, Ts,curr_N_sources);    % 1 - original, 2 - cyclic chirp, 3 - random power
        % Channel matrix A
        A = exp(1j*pi*[0:N_sensors-1].'*sin(theta));  
        X = A*s.';
        % AWGN
        sigma = 10^(-curr_SNR/20);
        w = sigma*(randn(N_sensors,length(n)) + 1j*(randn(N_sensors,length(n))))/sqrt(2); % noise

        X=X+w;        
        %% Original Fast-IF
        [IFF,ss] = Multi_Sensor_FAST_IF(X,N_sensors,65, N_sources, 6,100,0,0);
        
        y1 = 0;
        for iii=1:N_sources
            for jjj=1:N_sensors
                a(jjj,:)=ss(jjj,iii,:);
            end
            
            p=TMMUSIC(cov(a'), 2, N_sensors, 1, 1, theta_range');
            [x,y]=max(p);
            y1(iii)=y;
        end
        
        y1=y1-90;
        
        MSE_FAST_IF(curr_scenario)=mean((sort(y1/10)-sort(theta_deg/10)).^2);
        
        %% Adaptive Fast-IF
        optimal_th = 0.1;
        [IFF,ss,N_sources_est] = Adaptive_sources_FAST_IF(X,N_sensors,65, N_sensors, 6,100,optimal_th,0);
        
        y1 = 0;
        for iii=1:N_sources_est
            for jjj=1:N_sensors
                a(jjj,:)=ss(jjj,iii,:);
            end
            
            p=TMMUSIC(cov(a'), 2, N_sensors, 1, 1, theta_range');
            [x,y]=max(p);
            y1(iii)=y;
        end
        
        if(N_sources_est<N_sources)
            y1 = [y1,zeros(1,N_sources-N_sources_est)];
        elseif (N_sources_est>N_sources)
            theta_deg = [theta_deg, -90*ones(1,N_sources_est-N_sources)];
        end
            
        
        y1=y1-90;
        
        MSE_Adaptive_FAST_IF(curr_scenario)=mean((sort(y1/10)-sort(theta_deg/10)).^2);
    end
    MSE_FAST_IF_vec(SNR_idx)=mean(MSE_FAST_IF)
    MSE_Adaptive_FAST_IF_vec(SNR_idx)=mean(MSE_Adaptive_FAST_IF)
end

% FROM SIMULATION DONE in "Novel direction of arrival estimation using
% spatial adaptive"
plot(SNR,MSE_Adaptive_FAST_IF_vec,'--md','linewidth',2);
hold on;
plot(SNR,MSE_FAST_IF_vec,'r','linewidth',3);
hold off;
xlabel('Signal to Noise Ratio');
ylabel('Mean Square Error (dB)');
legend('Adaptive Fast IF','The FAST IF');
