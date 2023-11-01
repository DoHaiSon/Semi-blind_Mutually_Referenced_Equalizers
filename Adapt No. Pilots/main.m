clear all
close all
clc

num_exp   = 100;
Nt        = 2;         % number of tx antennas
Nr        = 4;         % number of rx antennas
M         = 3;         % channel order
N         = 10;        % window size
K         = M + N;     % rank of H
samp_size = 256;       % sample size
SNR       = -10:5:20;  % Signal to noise ratio (dB)
Nt_plot   = Nt;        % Number of tx in figure (max = Nt)

snr_i     = 12;
lambda    = 0.1;
N_p       = N;
N_p_f     = [N_p];
N_p_rc    = N;
N_p_rc_f  = [N_p_rc];
log_BER   = [];
log_BER_rc= [];
loss_BER  = [];
loss_BER_rc=[];
target_BER = 10^-4;
target_BER_rc = 10^-2;
iter = 1:100;
delta = 1;
delta_rc = 1;
trigger_1 = true;
trigger_2 = false;
trigger_rc = true;

%% Generate channel
% assumption i.i.d channels
% for tx = 1:Nt
%     H(tx, :, :) = (randn(Nr, M + 1) + 1i*randn(Nr, M + 1)) / sqrt(2);
% end
load H.mat

for i = 1:iter(end)

    fprintf('Iter %d \n', i);

    BER     = SB_MRE_adapt_func(N_p, H, num_exp, Nt, Nr, M, N, K, samp_size, Nt_plot, snr_i, lambda);
%     BER_rc  = SB_MRE_rc_adapt_func(N_p_rc, H, num_exp, Nt, Nr, M, N, K, samp_size, Nt_plot, snr_i, lambda);

    BER = mean(BER);
%     BER_rc = mean(BER_rc);

    if BER == 0
        BER = 10^-6;
    end

%     if BER_rc == 0
%         BER_rc = 10^-6;
%     end

    log_BER = [log_BER, BER];
%     log_BER_rc = [log_BER_rc, BER_rc];

    BER = log10(BER) - log10(target_BER);
%     BER_rc = log10(BER_rc) - log10(target_BER_rc) ;

    loss_BER = [loss_BER; BER];
%     loss_BER_rc = [loss_BER_rc; BER_rc];

    if loss_BER(end) == 0
        dx = 0;
    else
        dx = loss_BER(end);
    end
    
    % Adaptive 
    if ceil(dx) < 0
        if ~trigger_1 && trigger_2
            delta = 1;
            trigger_1 = true;
            trigger_2 = false;
        else
            delta = 2*delta;
        end
    else
        if trigger_1 && ~trigger_2
            delta = 1;
            trigger_1 = false;
            trigger_2 = true;
        else
            delta = 2*delta;
        end
    end
%     if ceil(dx) < 0
%         delta = 1;
%     else
%         delta = 2;
%     end
    N_p = N_p + delta * ceil(dx);

    if (N_p > samp_size)
        N_p = samp_size;
    end
    N_p_f    = [N_p_f; N_p];

%     if loss_BER_rc(end) == 0
%         dx = 0;
%     else
%         dx = loss_BER_rc(end);
%     end
% 
%     % Adaptive 
%     if ceil(dx) < 0
%         delta_rc = 1;
%     else
%         delta_rc = 2*delta;
%     end
%     N_p_rc = N_p_rc + delta_rc * ceil(dx);
% 
%     if (N_p_rc > samp_size)
%         N_p_rc = samp_size;
%     end
%     N_p_rc_f = [N_p_rc_f; N_p_rc];
end

N_p_f(end) = [];
N_p_rc_f(end) = [];
semilogy(iter, log_BER,    'LineWidth', 1);
hold on;
% semilogy(iter, log_BER_rc, 'LineWidth', 1);
% hold on;

% legend({'SB', 'SB reduced cost'});
legend({'SB'});
xlabel('Iteration');
ylabel('SER (dB)');
grid minor;