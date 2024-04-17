function [y,e] = varx_simulate(B,A,x,arg4)

% y = varx_simulate(B,A,x) Filters the vector input x through a VARX system
% with no noise. This is equivalent to a vectorial ARMA filter.
%
% y(t) = A*y(t-1) + B*x(t) 
% 
% here *y(t-1) and *x(t) represent concolution with the history of output
% an input with history na and nb duration, and A and B are matrizes of
% size [na,ydim,ydim] and [nb,ydim,xdim] respectively. Note that the
% current input x(t) does not affect the current output y(t). To implement
% instant effects, simply advance the input by one sample. x is the input
% of size [T,xdim] and y has size [T,ydim]. This option can be used to
% compute the impulse response of the VARX model, my calling it separatelly
% for each xdim set to an impulse.
%
% [y,e] = varx_simulate(B,A,x,se) if the third argument is a scalar, then
% computes output y(t) of a VARX system with external input x(t) and
% simulated innovation process e(t) with standard deviation se:
%
% y(t) = A*y(t-1) + B*x(t) + e(t)
%
% e(t) is size [T,ydim] and is drawn as zero mean normal i.i.d. noise.
%
% [yest,e] = varx_simulate(B,A,x,y) If the output y(t) is given, then the
% equation error is computed instead:
% 
% yest(t) = A*y(t-1) + B*x(t); 
% e(t) = yest(t) - y(t)
%
% If the intention is to simulate an AR model only, then set the input
% x=zeros(T,0); B=zeros(0,ydim,0); with xdim=0 and nb=0.
%
% See varx.m for how to estimate A and B from x,y data.

% (c) 08/01/2023 Lucas C Parra
%     08/02/2023 added zero padding at the start
%     08/15/2023 differing computation for equation error vs simulated noise 
%     08/16/2023 changed A, B format to be 3D tensor
%     08/24/2023 add option to compute ARMA filter putput with no innovation. Use tensorprod to avoid for loops
%     08/31/2023 changed input history to include current sample x(t) as predictor
%     09/13/2023 corrected the sign of the residual error to e=y-yest

% check all dimension
if size(A,2) ~= size(A,3) || size(A,2) ~= size(B,2) || size(B,3) ~= size(x,2)
    error('Dimentions of arguments inconsistent.')
end

[T,xdim] = size(x);
ydim = size(B,2);

% initialize with the additive innovation, or zero when computing equation
% error of just AMRA filtering
if ~exist('arg4','var') || isempty(arg4) % 4th argument missing: just filter
    e = zeros(T,ydim);
    y = e; 
    ER=0;
elseif size(arg4,1)==1 % 4th argument scalar: simulated additive innovation   
    e = arg4*randn(T,ydim); 
    y = e; 
    ER = 0; 
else % 4th argument vector: compute equation error
    %  this option made the code really uggly. Concider factoring out in
    %  the future as part of filterMIMO.m where it was before, and instead
    %  allow user to provided their own e(t) data here.
    y = zeros(T,ydim); 
    d = arg4; % called desired output when fitting B and A to observed x,y
    ER=1; 
end 

nb = size(B,1);
na = size(A,1);

% add zero padding
Q = max(na,nb);
x = [zeros(Q,xdim); x]; y = [zeros(Q,ydim); y]; if ER, d = [zeros(Q,ydim); d]; end

for t=Q+1:size(y,1)
    if ER % equation error based on observed output ycatual
        ypast = d(t-1:-1:t-na,:);
    else % simulate with random inovation and AR feedback
        ypast = y(t-1:-1:t-na,:);
    end
    xpast = x(t:-1:t-nb+1,:);
    y(t,:) = y(t,:) + tensorprod(A,ypast,[1 3],[1 2])' + tensorprod(B,xpast,[1 3],[1 2])';
end

% remove zero padding
y = y(Q+1:end,:); if ER, d = d(Q+1:end,:); end 

% either observed error, or simulated error
if ER, e = d - y; end