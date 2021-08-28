    function [theta_est,phi_est, N_Sources_est] = FAST_IF_2D_DOA(z_Sig,x_Sig, N_Sensors,N_Signals,win_length,delta,L,relative_energy_th,plot_flag,theta_phi_ref)
%% We present a 2D DOA estimation based on FAST-IF method:
% Inputs:
%   z_Sig - the signal sampled by the z-plain antennas
%   x_Sig - the signal sampled by the x-plain antennas
%   N_Sensors - number of sensors
%   N_Signals - number of sources
%   win_length - window length for FRFT
%   delta - the delta parameter for scanning frequencies
%   L - quantization parameter for scanning FRFT parameters
%   relative_energy_th - threshold for detection mode of algorithm
%   plot_flag - if 1 - plot estimations of algorithm
%   theta_phi_ref - vector of [theta,phi] of the real DOA
% 
% Outputs:
%   Fast_IF_Est_Vec - the instantaneous frequencies signals estimated
%   theta_est - Angles of arrival estimated in the Z plain
%   phi_est - Angles of arrival estimated in the X-Y plain
%   N_sources - Number of sources estimated (in detection mode)

%%
% for plot needs
ang_range = 0:1:180;
MUSIC_TH = 1.75;

% check if detection mode
if relative_energy_th == 0
    N_Sources_est = N_Signals;
    N_Signals_new = N_Signals;
else
    N_Sources_est = 0;
    N_Signals_new = N_Sensors;
end

%% Calculate angle of z plain - theta
[Fast_IF_est_z,Fast_IF_signals, N_sources_est] = Adaptive_sources_FAST_IF(z_Sig,N_Sensors,win_length, N_Signals_new, delta,L,relative_energy_th,0);

if plot_flag
    figure; 
    plot(Fast_IF_est_z'); title('Fast-IF Signal Estimation'); xlabel('Time [sec]'); ylabel('Frequency [Hz]');
end

for src_idx=1:N_Signals
    for curr_sensor=1:N_Sensors
        Signals_to_DOA(curr_sensor,:)=Fast_IF_signals(curr_sensor,src_idx,:);
    end
    
    p=TMMUSIC(cov(Signals_to_DOA'), 2, N_Sensors, 1, 1, (ang_range-90)');
    P(src_idx,:)=p/max(p);
    [max_peak, peak_idx] = max(p);
    N_sources_est_final = N_sources_est;
    % if found very close peak
    if max_peak < MUSIC_TH
        N_sources_est_final = src_idx-1;
%         MUSIC_TH - max_peak
        break;
    end
    [~,theta_est(src_idx)] = max(p);
    theta_est(src_idx) = 181 - theta_est(src_idx);
    theta_est_rad(src_idx) = theta_est(src_idx)*pi/180;
end
N_Signals;
% N_sources_est_final - N_Signals
% MUSIC_TH - max_peak
if plot_flag
    theta_ref = zeros(1,181);
    theta_ref(theta_phi_ref(1)+1)=1;
    figure;
    for src_idx = 1:size(P,1)
        plot(ang_range,fliplr(P(src_idx,:))); hold on; stem(ang_range,theta_ref); hold off;
    end 
    title('FAST-IF DOA Z-AXIS (\theta)'); xlabel('Degrees'); ylabel('Spatial Spectrum');
end

%% Calculate angle of x-y plain - phi
[Fast_IF_est_x,Fast_IF_signals] = Adaptive_sources_FAST_IF(x_Sig,N_Sensors,win_length, N_Signals, delta,L,relative_energy_th,0);

if plot_flag
    figure; 
    plot(Fast_IF_est_x'); title('Fast-IF Signal Estimation'); xlabel('Time [sec]'); ylabel('Frequency [Hz]');
end

for src_idx=1:N_Signals
    for curr_sensor=1:N_Sensors
        Signals_to_DOA(curr_sensor,:)=Fast_IF_signals(curr_sensor,src_idx,:);
    end
    
    p=TMMUSIC(cov(Signals_to_DOA'), 2, N_Sensors, 1, 1, (ang_range-90)');
    P(src_idx,:)=p/max(p);
    [~,gamma(src_idx)] = max(p);
    gamma(src_idx) = 181 - gamma(src_idx);
    gamma_rad(src_idx) = gamma(src_idx)*pi/180;
    cos_phi_est(src_idx) = 0.5 * (1 + cos(gamma_rad(src_idx)).^2-sin(gamma_rad(src_idx)).^2)...
                                    ./(cos(gamma_rad(src_idx)).*sin(theta_est_rad(src_idx)));  
    if cos_phi_est(src_idx) < -1
        cos_phi_est(src_idx) = cos_phi_est(src_idx) + pi;
    elseif cos_phi_est(src_idx) > 1
        cos_phi_est(src_idx) = cos_phi_est(src_idx) - pi;
    end
    phi_est_rad(src_idx) = acos(cos_phi_est(src_idx));
    phi_est(src_idx) = round(phi_est_rad(src_idx)*180/pi);
end

if plot_flag
    phi_ref = zeros(1,181);
    phi_ref(theta_phi_ref(2)+1)=1;
    figure;
    for src_idx = 1:size(P,1)
        plot(ang_range,fliplr(P(src_idx,:))); hold on; stem(ang_range,phi_ref); hold off;
    end
    title('FAST-IF DOA XY-PLAIN (\phi)'); xlabel('Degrees'); ylabel('Spatial Spectrum');
end

end


