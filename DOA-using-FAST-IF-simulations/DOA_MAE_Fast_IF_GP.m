% close all;
clear all;
addpath('C:\Program Files\MATLAB\R2021a\toolbox\tfsa_5_5');

%% Parameters:
N_sensors=4;
% samples of sig
Ts = 1;
fs = 1/Ts;
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

test_type = 2;

switch test_type
    case 1
        %% Regular SNR Test
        % Create signal sources
        [s, IF_O, N_sources] = Create_FM_Signals(1, n, Ts);    % 1 - original, 2 - cyclic chirp

        SNR_range = 20:2:30;   % dB
        test_length = 1000; % num of loops per SNR
        % Channel matrix A 
        A = exp(1j*pi*[0:N_sensors-1].'*sin(theta));                
        X = A*s.';                             % mixed source
        for SNR_idx = 1:length(SNR_range)    
            for curr_scenario=1:test_length 
                % AWGN     
                sigma = 10^(-SNR_range(SNR_idx)/20);
                w = sigma*(randn(N_sensors,length(n)) + 1j*(randn(N_sensors,length(n))))/sqrt(2); % noise       
                X=X+w;     

                [Fast_IF_IF,Fast_IF_signals] = Multi_Sensor_FAST_IF(X,N_sensors,65, N_sources, 4,100,0,0);

                for curr_src=1:N_sources
                    for curr_sensor=1:N_sensors
                        Signals_to_DOA(curr_sensor,:)=Fast_IF_signals(curr_sensor,curr_src,:);
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
            MAE_FAST_IF_vec(1,SNR_idx)=mean(MEAN_ABS_ERR_FAST_IF(1,:))
            MAE_FAST_IF_vec(2,SNR_idx)=mean(MEAN_ABS_ERR_FAST_IF(2,:))
        end

        figure;
        plot(SNR_range,MAE_FAST_IF_vec(1,:),'linewidth',2); grid; hold on;
        plot(SNR_range,MAE_FAST_IF_vec(2,:),'linewidth',2); hold off;
        title('Fast-IF Mean Absolute Error Vs SNR'); legend('1st Signal','2nd Signal'); 
        xlabel('SNR [dB]'); ylabel('Mean Absolute Error [degrees]');
        
    case 2
        %% Compare to adaptive
        % Create signal sources
        [s, IF_O, N_sources] = Create_FM_Signals(1, n, Ts);    % 1 - original, 2 - cyclic chirp

        SNR_range = 0:2:20;   % dB
        test_length = 1000; % num of loops per SNR
        % Channel matrix A 
        A = exp(1j*pi*[0:N_sensors-1].'*sin(theta));                
        X = A*s.';                             % mixed source
        for SNR_idx = 1:length(SNR_range)    
            curr_SNR = SNR_range(SNR_idx);
            for curr_scenario=1:test_length 
                % Generate signal     
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
                %% Fast IF
                [Fast_IF_IF,Fast_IF_signals] = Multi_Sensor_FAST_IF(X,N_sensors,65, N_sources, 6,100,0,0);
                y1 = 0;
                for curr_src=1:N_sources
                    for curr_sensor=1:N_sensors
                        Signals_to_DOA(curr_sensor,:)=Fast_IF_signals(curr_sensor,curr_src,:);
                    end

                    p=TMMUSIC(cov(Signals_to_DOA'), 2, N_sensors, 1, 1, theta_range');
                    [x,y]=max(p);
                    y1(curr_src)=y;
                end
                abs_err = abs((sort(y1)-sort(theta_deg)));
                MEAN_ABS_ERR_FAST_IF(1,curr_scenario)=abs_err(1);
                MEAN_ABS_ERR_FAST_IF(2,curr_scenario)=abs_err(2);
                
                %% Adaptive Fast IF
                [Fast_IF_IF,Fast_IF_signals] = Multi_Sensor_FAST_IF(X,N_sensors,65, N_sources, 6,100,0,0);
                y1 = 0;
                for curr_src=1:N_sources_est
                    for curr_sensor=1:N_sensors
                        Signals_to_DOA(curr_sensor,:)=Fast_IF_signals(curr_sensor,curr_src,:);
                    end

                    p=TMMUSIC(cov(Signals_to_DOA'), 2, N_sensors, 1, 1, theta_range');
                    [x,y]=max(p);
                    y1(curr_src)=y;
                end
                if(N_sources_est<N_sources)
                    y1 = [y1,zeros(1,N_sources-N_sources_est)];
                elseif (N_sources_est>N_sources)
                    theta_deg = [theta_deg, -90*ones(1,N_sources_est-N_sources)];
                end
                y1=y1-90;
                y1=y1-90;
                abs_err = abs((sort(y1)-sort(theta_deg)));
                MEAN_ABS_ERR_Adaptive_FAST_IF(1,curr_scenario)=abs_err(1);
                MEAN_ABS_ERR_Adaptive_FAST_IF(2,curr_scenario)=abs_err(2);

            end
            MAE_FAST_IF_vec(1,SNR_idx)=mean(MEAN_ABS_ERR_FAST_IF(1,:))
            MAE_FAST_IF_vec(2,SNR_idx)=mean(MEAN_ABS_ERR_FAST_IF(2,:))
            MAE_Adaptive_FAST_IF_vec(1,SNR_idx)=mean(MEAN_ABS_ERR_FAST_IF(1,:))
            MAE_Adaptive_FAST_IF_vec(2,SNR_idx)=mean(MEAN_ABS_ERR_FAST_IF(2,:))
        end
end

figure;
plot(SNR_range,MAE_FAST_IF_vec(1,:),'linewidth',2); grid; hold on;
plot(SNR_range,MAE_FAST_IF_vec(2,:),'linewidth',2); hold on;
plot(SNR_range,MAE_Adaptive_FAST_IF_vec(1,:),'linewidth',2); grid; hold on;
plot(SNR_range,MAE_Adaptive_FAST_IF_vec(2,:),'linewidth',2); hold on;
title('Fast-IF Mean Absolute Error Vs SNR'); 
legend('Regular Fast IF 1st Signal','Regular Fast IF 2nd Signal','Adaptive Fast IF 1st Signal','Adaptive Fast IF 2nd Signal'); 
xlabel('SNR [dB]'); ylabel('Mean Absolute Error [degrees]');
        
