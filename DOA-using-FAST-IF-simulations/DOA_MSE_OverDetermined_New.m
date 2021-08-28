close all;
clear all;
addpath('C:\Program Files\MATLAB\R2021a\toolbox\tfsa_5_5');

%% Parameters:
N_sources = 2;
N_sensors=3;
SNR=-10:2:10;
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
for i = 1:N_sources
    s_tmp(:,i) = exp(2*pi*1i*(init_freqs(i)*n+curve_slope(i)*n.^3/(N_samples*N_samples*3))).';
    s(:,i) = s_tmp(sel_samples,i);
    IF_O_tmp(i,:) = init_freqs(i)+curve_slope(i)*3*n.^2/(N_samples*N_samples*3);
    IF_O(i,sel_samples) = IF_O_tmp(i,sel_samples);
end

perc=0.4;

N_scenarios=100;
index=0;
delta=2;
for SNR_idx=1:length(SNR)
    curr_SNR = SNR(SNR_idx);
    % Channel matrix A 
    A = exp(1j*pi*[0:N_sensors-1].'*sin(theta));                
    X = A*s.';                             % mixed source
    for curr_scenario=1:N_scenarios 
        tic
        % AWGN     
        sigma = 10^(-curr_SNR/20);
        w = sigma*(randn(N_sensors,length(n)) + 1j*(randn(N_sensors,length(n))))/sqrt(2); % noise       
        X=X+w;
     
        [IF_tmp ss]= multi_sensor_source_separation_ridge_tracking_m(X, N_sources, delta,N_sensors);
        for iii=1:N_sources
            for jjj=1:N_sensors
                    a(jjj,:)=ss(jjj,iii,:);
            end
            
            p=TMMUSIC(cov(a'), 2, N_sensors, 1, 1, theta_range');
            [x,y]=max(p);
            y1(iii)=y;
        end
        
        y1=y1-90;
        MSE_RIDGE_TRACKING(curr_scenario)=mean((sort(y1/10)-sort(theta_deg/10)).^2);

        [IFF,ss] = Multi_Sensor_FAST_IF(X,N_sensors,65, N_sources, 2,100,0,0);

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
             
        %original code commented
%         tic
%         [ss,IF_out] = multi_sensor_source_separation_spatial_TF_direction(X, N_sources, 3,N_sensors);
%           %      [IFF,ss] = Multi_Sensor_FAST_IF(X,N_sensors,65, n_sources, 2,100,0,0);
%         toc
%         for iii=1:N_sources
%             for jjj=1:N_sensors
%                 a(jjj,:)=ss(jjj,iii,:);
%             end
%             
%             p=TMMUSIC(cov(a'), 2, N_sensors, 1, 1, theta_range');
%             [x,y]=max(p);
%             y1(iii)=y;
%         end
%         y1=y1-90;
%         mmssee_sadtfd(curr_scenario)=mean((sort(y1/10)-sort(theta_deg/10)).^2);
        MSE_SADTFD = 0;
        
        D   = mtfd(X, 'ckd',1, 0.3, 0.3, length(X));
        %%% Averaged Auto-TFD
        D_avg = zeros(length(X), length(X));
        for mm = 1:N_sensors, D_avg = D{mm,mm} + D_avg; end
        D_avg = D_avg./N_sensors;
        %%% Selection of high-energy (t,f) points
        thr = perc*max(max(D_avg));
        Tr = abs(D_avg) >= thr;
        [F_trace, ~] = find(Tr);
        n_p = length(F_trace);
        D_s = zeros(N_sensors, N_sensors);
        for m1 = 1:N_sensors
            for m2 = 1:N_sensors
                D_s(m1,m2) = (1/n_p).*sum(sum(D{m1,m2}.*Tr));
            end
        end

        %%% DOA Estimation
        P_tf_music_ckd = tf_music(D_s, N_sources, N_sensors, 2,1, theta_range);
       % P_tf_music_ckd =TMMUSIC(D_s, 2, N_sensors, n_sources, 1, theta_range);
        [~,y1]=findpeaks(P_tf_music_ckd,'NPeaks',4,'MinPeakDistance',10,'Threshold',0.001);
        y1=y1-90;
        %figure;plot(P_tf_music_ckd)
        if length(y1)<N_sources
            y1(length(y1)+1:N_sources)=0;
        elseif length(y1)>N_sources
            y1=y1(1:N_sources);
        end
      
        MSE_TF_MUSIC(curr_scenario)=mean((sort(y1/10)-sort(theta_deg/10)).^2);
%         toc
    end
    index=index+1;
    %mean(mmssee)
    MSE_FAST_IF_vec(index)=mean(MSE_FAST_IF)
    MSE_SADTFD_vec(index)=mean(MSE_SADTFD)
    MSE_TF_MUSIC(index)=mean(MSE_TF_MUSIC)
    MSE_RIDGE_TRACKING_vec(index)=mean(MSE_RIDGE_TRACKING)
end

% FROM SIMULATION DONE in "Novel direction of arrival estimation using
% spatial adaptive"
%snr_mse_sadtfd=[0.4418    0.0994    0.0615    0.0414    0.0272    0.0204];
plot(SNR,MSE_SADTFD_vec,'--md','linewidth',2);
hold on;
plot(SNR,MSE_FAST_IF_vec,'r','linewidth',3);
hold on;
plot(SNR,snr_mse_tf_music,'g','linewidth',3);
hold on;
%plot(SNR,10*(log10(snr_mse_post_proc)),'y','linewidth',2);
%hold on;
plot(SNR,MSE_RIDGE_TRACKING_vec,'b:','linewidth',3);

%openfig('DOA_MSE_over_determined -10 db to 10 db.fig')
xlabel('Signal to Noise Ratio');
ylabel('Mean Square Error (dB)');
legend('Adaptive Spatial TFDs','The FAST IF','Time-frequency Music','DOA based on IF estimation using ridge tracking');
%legend('The Proposed Method','Time-frequency Music','DOA based on IF estimation using ridge tracking');
