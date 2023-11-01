function BER = SB_MRE_adapt_func(N_p, H, num_exp, Nt, Nr, M, N, K, samp_size, Nt_plot, snr_i, lambda)

    err_sb    = zeros(1, Nt);
    BER       = [];

    for exp = 1:num_exp

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
        %% Blind MRE program
        X_n     = X(:, 1:samp_size-N);
        X_n_1   = X(:, 2:samp_size-N+1);
        A       = kron([eye((K-1)*Nt), zeros((K-1)*Nt, Nt)], X_n') - ...
            kron([zeros((K-1)*Nt, Nt), eye((K-1)*Nt)], X_n_1');
        Q       = A' * A;                                                   % N_rN KN_t x NrN KN_t

        %% Reduce the cost
        X_n     = X(:, 1: samp_size - N - K + 2);
        X_n_k   = X(:, K: samp_size - N + 1);
        A_rc    = kron([eye(Nt), zeros(Nt, Nt)], X_n') - ...
            kron([zeros(Nt, Nt), eye(Nt)], X_n_k');
        Q_rc    = A_rc' * A_rc; 

        %% ----------------------------------------------------------------
        %% Son's Semi-Blind program 
%         N_p     = N;
        X_SB    = X(:, 1:N_p - N + 1);                                      % N_tN x (N_p - N + 1)

        [G_sb, ~] = SB_MRE(X_SB, sig_src, Nt, Nr, N, K, M, N_p, Q, Q_rc, lambda);

        
        %% ----------------------------------------------------------------
        % Select the delay-th Equalizer
        delay       = round(K/2);

        % Equalization
        s_hat_sb    = G_sb   (:, :, delay)'    * X;
        s           = sig_src(:, K-delay+1:samp_size+M-delay+1);

        for tx  = 1:Nt_plot
           est_s_sb     = mydemod(s_hat_sb(tx, :));

           % Compare to src signals
           err_sb(tx)   = err_sb(tx)    + sum(s(tx, :)    ~= est_s_sb);
        end

        BER     = [BER; err_sb(tx)    / num_exp / samp_size];
    end
end