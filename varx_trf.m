function H = varx_trf(B,A,T)
% H = varx_trf(B,A,T) computes the MIMO finite impulse response of a varx
% model. In the context of neural responses to a stimulus a FIR is often
% refered to as "temporal response function" and the MIMO version is the
% multivariate RRF (mTRF). T is the desired length of the FIR. ARMA filters
% B and A are organized as in the varx() function.

% (c) August 30, 2023 Lucas C Parra

xdim = size(B,3);
for i=1:xdim
    x = zeros(T,xdim);
    x(1,i) = 1;
    H(:,:,i) = varx_simulate(B,A,x);
end
