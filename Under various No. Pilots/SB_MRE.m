function [G, G_rc] = SB_MRE(X_SB, S_n, Nt, Nr, N, K, M, N_p, Q, Q_rc)
    %% Declear variables
    lambda  = 0.1;  

    S_SB    = [];                                             % N_tK x (N_p - N + 1)
    S_rc    = [];                                             % N_tP x (N_p - N + 1)
    for i   = N:N_p                                           % S = [S_{1}(N-1), ..., S_{1}(N_p-1); S_{N_t}(N-1), ..., S_{N_t}(N_p-1)]
        S_k = S_n(:, i+M:-1:i+M-K+1);
        S_SB= [S_SB, S_k(:)];

        S_rc_k = [S_n(:, i+M) S_n(:, i+M-K+1)];
        S_rc   = [S_rc, S_rc_k(:)];
    end

    %% J(G)                         
    A       = kron(eye(K*Nt), X_SB');                         % KN_t(N_p - N + 1) x KN_tN_rN
    s_bar   = S_SB';                                          
    s_bar   = s_bar(:);                                       % KN_t(N_p - N + 1) x 1
    
    g       = pinv(A' * A + lambda * Q) * A' * s_bar;         % KN_rNN_t x 1
    % Reshape G
    G       = reshape(g, [Nr*N, Nt, K]);                      % N_rN x N_t x K

    %% J(G_rc)
    A       = kron(eye(2*Nt), X_SB');                         % 2N_t(N_p - N + 1) x 2N_tN_rN
    s_bar   = S_rc';                                          
    s_bar   = s_bar(:);                                       % 2N_t(N_p - N + 1) x 1
    
    g_rc    = pinv(A' * A + lambda * Q_rc) * A' * s_bar;      % 2N_rNN_t x 1
    % Reshape G
    G_rc    = reshape(g_rc, [Nr*N, Nt, 2]);                   % N_rN x N_t x 2
end