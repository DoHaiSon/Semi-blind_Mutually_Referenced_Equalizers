function [est] = mydemod(s_hat)
    est = ( sign(real(s_hat)) + 1i*sign(imag(s_hat)) ) / sqrt(2);
end