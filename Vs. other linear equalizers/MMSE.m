function s_hat = MMSE(X, H, snr)    
    eye_shape = size(H);
%     sigma_2 = mean(mean(noise * noise'));              % noise power;
%     E_x     = mean(mean(X * X'));                      % avg power / symbol;

    s_hat   = inv(H' * H + (10^(-snr/10) * eye(eye_shape(2)))) * H' * X;
end