clear all
close all
clc

num_exp   = 50000;
Nt        = 2;         % number of tx antennas
Nr        = 4;         % number of rx antennas
M         = 3;         % channel order
N         = 10;        % window size
K         = M + N;     % rank of H
samp_size = 256;       % sample size
SNR       = -10:2:20;  % Signal to noise ratio (dB)
Nt_plot   = Nt;        % Number of tx in figure (max = Nt)

res         = [];
res_sb      = [];
res_rc      = [];
res_sb_rc   = [];
res_zf      = [];
res_mmse    = [];
res_mle     = [];

for snr_i = SNR
    fprintf('-------------------------------------------------------------\nWorking at SNR: %d dB\n', snr_i);
    err       = zeros(1, Nt);
    err_sb    = zeros(1, Nt);
    err_rc    = zeros(1, Nt);
    err_sb_rc = zeros(1, Nt);
    err_zf    = zeros(1, Nt);
    err_mmse  = zeros(1, Nt);
    err_mle   = zeros(1, Nt);

    for exp = 1:num_exp
        fprintf('Experience No %d \n', exp);
        %% Generate channel
        % assumption i.i.d channels
        for tx = 1:Nt
            H(tx, :, :) = (randn(Nr, M + 1) + 1i*randn(Nr, M + 1)) / sqrt(2);
        end
        
        H_toep = [];
        for tx = 1:Nt
            h_toep = [];
            Ch = squeeze(H(tx, :, :));
            for rx = 1:Nr
                padding      = zeros(1, samp_size);
                padding(1,1) = Ch(rx, 1);
                h            = [Ch(rx, :), zeros(1, samp_size - 1)];
                h_toep       = [h_toep; toeplitz(padding, h)];
            end
            H_toep           = cat(2, H_toep, h_toep);
        end

        %% Generate signals
        % Signal src (QPSK)
        for tx = 1:Nt
            sig_src(tx, :) = (sign(rand(samp_size + M,1)-0.5) + 1i*sign(rand(samp_size + M,1)-0.5))/sqrt(2);
        end
            
        S_shape= Nt * (samp_size + M);
        X_zf   = H_toep * reshape(sig_src, [S_shape, 1]);
        noise  = X_zf;
        X_zf   = awgn(X_zf, snr_i);
        noise  = X_zf - noise;

        % Signal rec
        sig_rec = zeros(samp_size + 2*M, Nr);
        for tx = 1:Nt
            y = [];
            for rx = 1:Nr
                y(:, rx) = conv( squeeze(H(tx, rx, :)).', sig_src(tx, :) ) ;
            end
            sig_rec = sig_rec + y;
        end

        sig_rec = sig_rec(M+1:samp_size + M, :);                               % samp_size x Nr

        sig_rec = awgn(sig_rec, snr_i);

        %% MRE
        X       = [];
        for ii  = 1:Nr
            x     = sig_rec(:, ii);
            mat   = hankel(x(1:N), x(N:samp_size));
            mat   = mat(N:-1:1, :);
            X     = [X; mat];
        end

        %% ----------------------------------------------------------------
        %% Non-Blind MRE program
        s_hat_zf = reshape(ZF(X_zf, H_toep), [Nt, samp_size + M]);
        s_hat_mmse = reshape(MMSE(X_zf, H_toep, snr_i), [Nt, samp_size + M]);
%         s_hat_mle = reshape(MLE(X_zf, H_toep, reshape(sig_src, [S_shape, 1])), [Nt, samp_size + M]);
        
        %% ----------------------------------------------------------------
        %% Blind MRE program
        X_n     = X(:, 1:samp_size-N);
        X_n_1   = X(:, 2:samp_size-N+1);
        A       = kron([eye((K-1)*Nt), zeros((K-1)*Nt, Nt)], X_n') - ...
            kron([zeros((K-1)*Nt, Nt), eye((K-1)*Nt)], X_n_1');
        Q       = A' * A;                                                   % N_rN KN_t x NrN KN_t
        [V_b, ~]= eigs(Q, 1,'SM');                                          % N_rN KN_t x 1

        % Reshape V full equalizers matrix
        V = reshape(V_b, [Nr*N, Nt, K]);     % short code

        %% Reduce the cost
        X_n     = X(:, 1: samp_size - N - K + 2);
        X_n_k   = X(:, K: samp_size - N + 1);
        A_rc    = kron([eye(Nt), zeros(Nt, Nt)], X_n') - ...
            kron([zeros(Nt, Nt), eye(Nt)], X_n_k');
        Q_rc    = A_rc' * A_rc; 
        [V_b_rc, ~]= eigs(Q_rc, 1,'SM');       

        % Reshape V full equalizers matrix
        V_rc    = reshape(V_b_rc, [Nr*N, Nt, 2]);     % short code

        %% ----------------------------------------------------------------
        %% Son's Semi-Blind program 
        N_p     = samp_size / 8;
        X_SB    = X(:, 1:N_p - N + 1);                                      % N_tN x (N_p - N + 1)

        [G_sb, G_sb_rc] = SB_MRE(X_SB, sig_src, Nt, Nr, N, K, M, N_p, Q, Q_rc, 0.1);

        
        %% ----------------------------------------------------------------
        % Select the delay-th Equalizer
        delay       = round(K/2);
        delay_rc    = 2;

        % Equalization
        s_hat       = V      (:, :, delay)'    * X;
        s_hat_rc    = V_rc   (:, :, delay_rc)' * X;
        s_hat_sb    = G_sb   (:, :, delay)'    * X;
        s_hat_sb_rc = G_sb_rc(:, :, delay_rc)' * X;
    
        s           = sig_src(:, K-delay+1:samp_size+M-delay+1);
        idelay      = delay_rc * (K - 1) - K + 2;   % 0 or (K-1)-th delay
        s_rc        = sig_src(:, K-idelay+1:samp_size+M-idelay+1);

        % Remove the inherent scalar indeterminacy
        s_hat       = s    * pinv(s_hat)     * s_hat;
        s_hat_rc    = s_rc * pinv(s_hat_rc)  * s_hat_rc;

        for tx  = 1:Nt_plot
           est_s        = mydemod(s_hat(tx, :));
           est_s_sb     = mydemod(s_hat_sb(tx, :));
           est_s_rc     = mydemod(s_hat_rc(tx, :));
           est_s_sb_rc  = mydemod(s_hat_sb_rc(tx, :));
           est_s_zf     = mydemod(s_hat_zf(tx, :));
           est_s_mmse   = mydemod(s_hat_mmse(tx, :));
%            est_s_mle    = mydemod(s_hat_mle(tx, :));

           % Compare to src signals
           err(tx)      = err(tx)       + sum(s(tx, :)    ~= est_s);
           err_sb(tx)   = err_sb(tx)    + sum(s(tx, :)    ~= est_s_sb);
           err_rc(tx)   = err_rc(tx)    + sum(s_rc(tx, :) ~= est_s_rc);
           err_sb_rc(tx)= err_sb_rc(tx) + sum(s_rc(tx, :) ~= est_s_sb_rc);
           err_zf(tx)   = err_zf(tx)    + sum(sig_src(tx, :) ~= est_s_zf);
           err_mmse(tx) = err_mmse(tx)  + sum(sig_src(tx, :) ~= est_s_mmse);
%            err_mle(tx)  = err_mle(tx)   + sum(sig_src(tx, :) ~= est_s_mle);
        end
    end
    
    for tx = 1:Nt_plot
        res      (tx, find(SNR == snr_i)) = err(tx)       / num_exp / samp_size;
        res_sb   (tx, find(SNR == snr_i)) = err_sb(tx)    / num_exp / samp_size;
        res_rc   (tx, find(SNR == snr_i)) = err_rc(tx)    / num_exp / samp_size;
        res_sb_rc(tx, find(SNR == snr_i)) = err_sb_rc(tx) / num_exp / samp_size;
        res_zf   (tx, find(SNR == snr_i)) = err_zf(tx)    / num_exp / samp_size;
        res_mmse (tx, find(SNR == snr_i)) = err_mmse(tx)  / num_exp / samp_size;
%         res_mle  (tx, find(SNR == snr_i)) = err_mle(tx)   / num_exp / samp_size;
    end
end

res = mean(res);
res_sb = mean(res_sb);
res_rc = mean(res_rc);
res_sb_rc = mean(res_sb_rc);
res_zf = mean(res_zf);
res_mmse = mean(res_mmse);

semilogy(SNR, res_zf,    'LineWidth', 1);
hold on;
semilogy(SNR, res_mmse,  'LineWidth', 1);
hold on;
semilogy(SNR, res,       'LineWidth', 1);
hold on;
semilogy(SNR, res_rc,    'LineWidth', 1);
hold on;
semilogy(SNR, res_sb,    'LineWidth', 1);
hold on;
semilogy(SNR, res_sb_rc, 'LineWidth', 1);
hold on;
%     semilogy(SNR, res_mle(tx, :), 'LineWidth', 1);

legend({'ZF', 'MMSE', 'B', 'B reduced cost', 'SB', 'SB reduced cost'});
xlabel('SNR (dB)');
ylabel('SER (dB)');
grid minor;