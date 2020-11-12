function [lpf] = lpf_FIR(omegaC, n)
    shifted_0 = (n-1)/2;
    M = (n-1);
    shifted_time = [0:1:M];
    lpf = sin(omegaC*(shifted_time-shifted_0+eps))./(pi*(shifted_time-shifted_0+eps));
end