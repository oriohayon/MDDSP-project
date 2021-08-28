% close all;
clear all;
% Need to have to following tool box:
addpath('C:\Program Files\MATLAB\R2021a\toolbox\tfsa_5_5');

%% Parameters:
N_sensors=4;
% samples of sig
Ts = 1;
fs = 1/Ts;
time_axis = 0:Ts:128-Ts;
N_samples = length(time_axis);
n=0:N_samples-1;
nT = n*Ts;
% generate noise
SNR=50;
theta_range = -90:1:90;

%% Randomize signal sources
for i = 2
switch (i)
    case 1
        %% Find optimal Threshold of energy
        N_scenarios = 1000;
        TH_range = 0.05:0.05:0.3;
        % TH_range = 100*TH_range;
        for TH_idx = 1:length(TH_range)
            curr_TH = TH_range(TH_idx);
            err_cnt = 0;
            err_cnt_1 = 0;
            err_cnt_2 = 0;   
            for curr_scenario = 1:N_scenarios
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

                % Additive WGN
                sigma = 10^(-SNR/20);
                w = sigma*(randn(N_sensors,length(n)) + 1j*(randn(N_sensors,length(n))))/sqrt(2); % noise
                X=X+w;

                %% Fast IF algorithm:
                [Fast_IF_IF,Fast_IF_signals,N_sources_est] = Adaptive_sources_FAST_IF(X,N_sensors,65, N_sensors, 6,100,curr_TH,curr_N_sources);
                if N_sources_est > curr_N_sources
                    err_cnt_1 = err_cnt_1 + 1;
                    err_cnt = err_cnt + 1;
                elseif N_sources_est < curr_N_sources
                    err_cnt_2 = err_cnt_2 + 1;
                    err_cnt = err_cnt + 1;
                end
            end
            Over_Estimate_Rate(TH_idx) = err_cnt_1/N_scenarios
            Under_Estimate_Rate(TH_idx) = err_cnt_2/N_scenarios
            Scenarios_Error_Rate(TH_idx) = err_cnt/N_scenarios
        end
        figure;
        plot(TH_range,Scenarios_Error_Rate);

        figure;
        plot(TH_range,Over_Estimate_Rate); title('Over Estimate Rate');

        figure;
        plot(TH_range,Under_Estimate_Rate); title('Under Estimate Rate');
        
%% Vs SNR
case 2 
        Over_Estimate_Rate = 0;
        Under_Estimate_Rate = 0;
        Scenarios_Error_Rate = 0;
        N_scenarios = 1000;
        N_batches_to_det = 4;
%         SNR_range = 4; 
        SNR_range = 0:2:24;
        curr_TH = 0.1;
        MUSIC_TH = 1.75;
        for SNR_idx = 1:length(SNR_range)
            SNR = SNR_range(SNR_idx);
            err_cnt = 0;
            err_cnt_1 = 0;
            err_cnt_2 = 0;   
            err_dist = 0;
                        
            for curr_scenario = 1:N_scenarios
                curr_N_sources = randi(N_sensors);
                mean_n_sources_est = 0;
                % randomize angle with +-4 degrees deviation
                theta = (randi(171,[1,curr_N_sources])-91)*pi/180;
                theta_deg=round(theta *180/pi) + (randi([-4,4],size(theta)));
                for curr_batch = 1:N_batches_to_det                  
                    [s, IF_O, N_sources] = Create_FM_Signals(3, n, Ts,curr_N_sources);    % 1 - original, 2 - cyclic chirp, 3 - random power
                    % Channel matrix A
                    A = exp(1j*pi*[0:N_sensors-1].'*sin(theta));  
                    X = A*s.';

                    % Additive WGN
                    sigma = 10^(-SNR/20);
                    w = sigma*(randn(N_sensors,length(n)) + 1j*(randn(N_sensors,length(n))))/sqrt(2); % noise
                    X=X+w;

                    %% Fast IF algorithm:
                    [Fast_IF_IF,Fast_IF_signals,N_sources_est, confidence_lvl] = Adaptive_sources_FAST_IF(X,N_sensors,65, N_sensors, 6,100,curr_TH,curr_N_sources);
                    N_sources_est_final = N_sources_est;
                    y1 = 0;
                    P = 0;
                    last_peaks = [];
                    for iii=1:N_sources_est
                        for jjj=1:N_sensors
                            a(jjj,:)=Fast_IF_signals(jjj,iii,:);
                        end

                        p=TMMUSIC(cov(a'), 2, N_sensors, 1, 1, theta_range');
                        [max_peak, peak_idx] = max(p);
                        if iii == 1
                            max_peak_ref = max_peak;
                        end
                        % if found low peak - remove estimation
                        if max_peak < MUSIC_TH
                            N_sources_est_final = N_sources_est_final - 1;
                        end

                        y1(iii)=peak_idx;
                    end
                    mean_n_sources_est = mean_n_sources_est + N_sources_est_final;
                end
                mean_n_sources_est = round(mean_n_sources_est/N_batches_to_det);
                N_sources_est_final = mean_n_sources_est;
                if N_sources_est_final > curr_N_sources
                    err_cnt_1 = err_cnt_1 + 1;
                    err_cnt = err_cnt + 1;
                    err_dist = err_dist + N_sources_est_final-curr_N_sources;
                    dummy = 0;
                elseif N_sources_est_final < curr_N_sources
                    err_cnt_2 = err_cnt_2 + 1;
                    err_cnt = err_cnt + 1;
                    err_dist = err_dist + curr_N_sources-N_sources_est_final;
                    dummy = 0;   
                end               
            end
                       
            Over_Estimate_Rate(SNR_idx) = err_cnt_1/N_scenarios
            Under_Estimate_Rate(SNR_idx) = err_cnt_2/N_scenarios
            Scenarios_Error_Rate(SNR_idx) = err_cnt/N_scenarios
            mean_err_dist(SNR_idx) = err_dist/err_cnt
        end

        figure;
        plot(SNR_range,Scenarios_Error_Rate); title('Total Error Rate');
        xlabel('SNR[dB]'); 
        figure;
        plot(SNR_range,Over_Estimate_Rate); title('Over Estimate Rate');
        xlabel('SNR[dB]');
        figure;
        plot(SNR_range,Under_Estimate_Rate); title('Under Estimate Rate');
        xlabel('SNR[dB]');
        figure;
        plot(SNR_range,mean_err_dist); title('Mean Absolute Error Distance');
        xlabel('SNR[dB]');
       
%% Vs MUSIC TH
case 3
    Over_Estimate_Rate = 0;
    Under_Estimate_Rate = 0;
    Scenarios_Error_Rate = 0;
    N_scenarios = 500;
    SNR = 6;
    curr_TH = 0.1;
%     MUSIC_TH_range = 1:0.25:2.5;
    MUSIC_TH_range = 0.01:0.01:0.1;
    for MUSIC_TH_idx = 1:length(MUSIC_TH_range)
        MUSIC_TH = MUSIC_TH_range(MUSIC_TH_idx);
        err_cnt = 0;
        err_cnt_1 = 0;
        err_cnt_2 = 0;   
        for curr_scenario = 1:N_scenarios
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

            % Additive WGN
            sigma = 10^(-SNR/20);
            w = sigma*(randn(N_sensors,length(n)) + 1j*(randn(N_sensors,length(n))))/sqrt(2); % noise
            X=X+w;

            %% Fast IF algorithm:
            [Fast_IF_IF,Fast_IF_signals,N_sources_est, sources_energy] = Adaptive_sources_FAST_IF(X,N_sensors,65, N_sensors, 6,100,curr_TH,0);
            N_sources_est_final = N_sources_est;
            y1 = 0;
            P = 0;
            last_peaks = [];
            for iii=1:N_sources_est
                for jjj=1:N_sensors
                    a(jjj,:)=Fast_IF_signals(jjj,iii,:);
                end

                p=TMMUSIC(cov(a'), 2, N_sensors, 1, 1, theta_range');
                [max_peak, peak_idx] = max(p);
                if iii == 1
                    max_peak_ref = max_peak;
                end
                % if found very close peak
                if max_peak < MUSIC_TH
                    N_sources_est_final = iii-1;
                    break;
                end
                y1(iii)=peak_idx;
            end
            if N_sources_est_final > curr_N_sources
                err_cnt_1 = err_cnt_1 + 1;
                err_cnt = err_cnt + 1;
            elseif N_sources_est_final < curr_N_sources
                err_cnt_2 = err_cnt_2 + 1;
                err_cnt = err_cnt + 1;
            end
        end
        Over_Estimate_Rate(MUSIC_TH_idx) = err_cnt_1/N_scenarios
        Under_Estimate_Rate(MUSIC_TH_idx) = err_cnt_2/N_scenarios
        Scenarios_Error_Rate(MUSIC_TH_idx) = err_cnt/N_scenarios
    end
    figure;
    plot(MUSIC_TH_range,Scenarios_Error_Rate); title('Total Error Rate');
    xlabel('MUSIC TH'); 
    figure;
    plot(MUSIC_TH_range,Over_Estimate_Rate); title('Over Estimate Rate');
    xlabel('MUSIC TH'); 
    figure;
    plot(MUSIC_TH_range,Under_Estimate_Rate); title('Under Estimate Rate');
    xlabel('MUSIC TH'); 
end
end
% figure; 
% plot(Fast_IF_IF'); hold on; plot(IF_O','color','r','LineStyle','-.','LineWidth',0.75);
% title('Fast-IF Signal Estimation'); xlabel('Time [sec]'); ylabel('Frequency [Hz]');
% figure;
% plot(theta_range,P')
% hold on; stem(theta_range,theta_axis)
% title('FAST-IF DOA'); xlabel('Degrees'); ylabel('Spatial Spectrum');
