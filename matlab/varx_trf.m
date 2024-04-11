function [H,Ainv] = varx_trf(B,A,T)
% H = varx_trf(B,A,T) computes the MIMO impulse response of a varx model.
% In the context of neural responses to a stimulus a impulse response is
% often refered to as "temporal response function" and the MIMO version is
% the multivariate RRF (mTRF). T is the desired length of the impulse
% response. ARMA filters B and A are organized as in the varx() function.
%
% [~,Ainv] = varx_trf(B,A,T) computes in Ainv the impulse response from the
% innovation process to the output of the VARX model. Essentially, how does
% a unit injection of innovation into recurrent dynamic affect the output
% variables. This corresponds to the inverse of A. If the VARX model has
% not input, i.e. its only a VAR model, then set B=zeros(T,ydim,0). 

% (c) August 30, 2023 Lucas C Parra
%     February 29, 2024 Added computation of Ainv

[~,ydim,xdim] = size(B);
H = []; % empty unless xdim>0
for i=1:xdim
    x = zeros(T,xdim);
    x(1,i) = 1;
    H(:,:,i) = varx_simulate(B,A,x);
end

if nargout>1
    for i=1:ydim
        x = zeros(T,1); x(1) = 1; 
        B = zeros(1,ydim,1); B(1,i) = 1;
        Ainv(:,:,i) = varx_simulate(B,A,x);
    end
end