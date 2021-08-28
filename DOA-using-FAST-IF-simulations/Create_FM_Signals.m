function [s, IF_O, N_sources] = Create_FM_Signals(signal_idx, n, Ts, N_sources)

nT = n*Ts;
if nargin < 4
    N_sources = 2;
end

switch signal_idx
    %% original chirp signal
    case 1  
        for i = 1:N_sources
            % Parameters for signals:
            power_norm = [1 , 1];
            init_freqs = [0.4 , 0.25];
            curve_slope = [-0.3 , 0.2]./(3*128*128);
            chirp_func_integral = curve_slope(i)*nT.^3;
            chirp_func = curve_slope(i)*3*nT.^2;
            % create FM signals
            s(:,i) = power_norm(i).*exp(2*pi*1i*(init_freqs(i)*nT+chirp_func_integral)).';
            IF_O(i,:) = init_freqs(i)+chirp_func;
        end
    %% Cyclic chirp    
    case 2
        % Parameters for signals:
        power_norm = [0.6 , 1];
        init_freqs = [0.1 , 0.4];
        slopes = [0.3, -3]./(3*128*128);
        % first signal
        chirp_func_integral = slopes(1)*nT.^3;
        chirp_func = slopes(1)*3*nT.^2;
        s(:,1) = power_norm(1).*exp(2*pi*1i*(init_freqs(1)*nT+chirp_func_integral)).';
        IF_O(1,:) = init_freqs(1)+chirp_func;
        % Second signal (cyclic chirp)
        n_cycles = 4;
        cycle_length = floor(length(nT)/n_cycles);
        chirp_func_integral_one_cycle = slopes(2)*nT(1:cycle_length).^3;
        chirp_func_integral = repmat(chirp_func_integral_one_cycle,1,n_cycles);
        chirp_func_integral = [chirp_func_integral,zeros(1,length(nT)-length(chirp_func_integral))];
        chirp_func_one_cycle = slopes(2)*3*nT(1:cycle_length).^2;
        chirp_func = repmat(chirp_func_one_cycle,1,n_cycles);
        chirp_func = [chirp_func,zeros(1,length(nT)-length(chirp_func))];
        s(:,2) = power_norm(2).*exp(2*pi*1i*(init_freqs(2)*nT+chirp_func_integral)).';
        IF_O(2,:) = init_freqs(2)+chirp_func;
    
    %% Load from stack and Randomize power of signals  
    case 3
        % Power randomize for signals:
%         power_norm = [1, 0.7.*rand(1,N_sources-1)+0.3];
        power_norm = ones(1,N_sources);
        load('IF_Stack');
        load('Signals_Stack');
        for i = 1:N_sources
            curr_idx = randi(N_sources);
            s(:,i) = power_norm(i).*Signals_Stack(curr_idx,:);
            Signals_Stack = [Signals_Stack(1:curr_idx-1,:);Signals_Stack(curr_idx+1:end,:)];
            IF_O(i,:) = IF_Stack(curr_idx,:);
            IF_Stack = [IF_Stack(1:curr_idx-1,:);IF_Stack(curr_idx+1:end,:)];
        end
    %% steps chirp signal
    case 4  
        % create steps functions:
        step_length = 32;
        last_integral = 0;
        for k = 0:step_length:N_samples-1
            step_func_1(k+1:k+step_length) = 0.5 - k/(2*N_samples);
            step_func_integral_1(k+1:k+step_length) = step_func_1(k+1)*n(1:step_length)+last_integral;
            last_integral = step_func_integral_1(k+step_length) + step_func_1(k+1);
        end

        last_integral = 0;
        for k = 0:step_length:N_samples-1
            step_func_2(k+1:k+step_length) = k/(2*N_samples);
            step_func_integral_2(k+1:k+step_length) = step_func_2(k+1)*n(1:step_length)+last_integral;
            last_integral = step_func_integral_2(k+step_length) + step_func_2(k+1);
        end

        steps_function = [step_func_1;step_func_2]/3;
        steps_function_integral = [step_func_integral_1;step_func_integral_2]/3;
        

end