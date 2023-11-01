function s_hat = ZF(X, H)    
    s_hat = inv(H' *H) * H' * X;
end