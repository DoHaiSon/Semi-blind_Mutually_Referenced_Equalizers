clear all
close all
clc

num_exp   = 2000;
Nt        = 2;         % number of tx antennas
Nr        = 4;         % number of rx antennas
M         = 3;         % channel order
N         = 10;        % window size
K         = M + N;     % rank of H
samp_size = 256;       % sample size
SNR       = 10;        % Signal to noise ratio (dB)

res_sb      = [];
res_sb_rc   = [];

P       = N:5:256;

for N_p = P
    fprintf('-------------------------------------------------------------\nWorking at N_p: %d pilots\n', N_p);
    err_sb    = zeros(1, Nt);
    err_sb_rc = zeros(1, Nt);

    for exp = 1:num_exp
        fprintf('Experience No %d \n', exp);

        %% Generate channel
        % assumption i.i.d channels
        H   = [];
        for tx = 1:Nt
            H(tx, :, :) = (randn(Nr, M + 1) + 1i*randn(Nr, M + 1)) / sqrt(2);
        end

        %% Generate signals
        % Signal src (QPSK)
        for tx = 1:Nt
            sig_src(tx, :) = (sign(rand(samp_size + M,1)-0.5) + 1i*sign(rand(samp_size + M,1)-0.5))/sqrt(2);
        end

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

        sig_rec = awgn(sig_rec, SNR);

        %% MRE
        X       = [];
        for ii  = 1:Nr
            x     = sig_rec(:, ii);
            mat   = hankel(x(1:N), x(N:samp_size));
            mat   = mat(N:-1:1, :);
            X     = [X; mat];
        end

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
        X_SB    = X(:, 1:N_p - N + 1);                                      % N_tN x (N_p - N + 1)

        [G_sb, G_sb_rc] = SB_MRE(X_SB, sig_src, Nt, Nr, N, K, M, N_p, Q, Q_rc);

        
        %% ----------------------------------------------------------------
        % Select the delay-th Equalizer
        delay       = round(K/2);
        delay_rc    = 1;

        % Equalization
        s_hat_sb    = G_sb   (:, :, delay)'    * X;
        s_hat_sb_rc = G_sb_rc(:, :, delay_rc)' * X;
    
        s           = sig_src(:, K-delay+1:samp_size+M-delay+1);
        idelay      = delay_rc * (K - 1) - K + 2;   % 0 or (K-1)-th delay
        s_rc        = sig_src(:, K-idelay+1:samp_size+M-idelay+1);

        for tx  = 1:Nt
           est_s_sb     = ( sign(real(s_hat_sb(tx, :)))     + 1i*sign(imag(s_hat_sb(tx, :))) )    / sqrt(2);
           est_s_sb_rc  = ( sign(real(s_hat_sb_rc(tx, :)))  + 1i*sign(imag(s_hat_sb_rc(tx, :))) ) / sqrt(2);

           % Compare to src signals
           err_sb(tx)   = err_sb(tx)    + sum(s(tx, :)    ~= est_s_sb);
           err_sb_rc(tx)= err_sb_rc(tx) + sum(s_rc(tx, :) ~= est_s_sb_rc);
        end
    end
    
    for tx = 1:Nt
        res_sb   (tx, find(P == N_p)) = err_sb(tx)    / num_exp / samp_size;
        res_sb_rc(tx, find(P == N_p)) = err_sb_rc(tx) / num_exp / samp_size;
    end
end

res_sb = mean(res_sb);
res_sb_rc = mean(res_sb_rc);

semilogy(P, res_sb,    'LineWidth', 1);
hold on;
semilogy(P, res_sb_rc, 'LineWidth', 1);

legend({'SB', 'SB reduced cost'});
xlabel('Number of pilot symbols');
ylabel('SER (dB)');
grid minor;