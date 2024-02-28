clear all
T=100; % try T not much bigger than xdim
n=[2 40 10]; % number of dimensions to null and test corresponding pvalues
xdim = sum(n);
gamma = 0.5/(T-xdim); 

% random simulations to estimate empirical pvalue 
for i=1000:-1:1
    x = randn(T,xdim);
    h = randn(xdim,1);
    h(1:n(1)) = 0; % null the first n to test pvalues
    h(sum(n(1:2))+1:end) = 0; % null the last n to test pvalues
    e = randn(T,1);
    y = x*h + e;
    pval(i,:) = simulate(x,y,gamma,n);
end
mean(pval<0.05) % should return [0.05 1.0 0.05] if pvalues are correctly estimated

function [pval] = simulate(x,y,gamma,n)

[s2,Bias] = fit(x,y,gamma);

[T,xdim] = size(x);
for i=length(n):-1:1 
    ii = [1:sum(n(1:i-1)) sum(n(1:i))+1:xdim]; % use these inputs excluding n inputs
    [s2r,Biasr] = fit(x(:,ii),y,gamma);
    D = (T-sum(n)+n(i))*log(s2r/s2) - Biasr + Bias; % not the exact formula, but I calibrated and seems to work well for small T
    pval(i) = 1-chi2cdf(D,n(i));
end
end

function [s2,Bias] = fit(x,y,gamma)

[T,xdim] = size(x);
Rxx = x'*x;
Rxy = x'*y;

% Ridge regularizer
I = eye(size(Rxx))*mean(real(eig(Rxx)));

% Least squares estimate with regularization
h = ((1-gamma)*Rxx+gamma*I)\Rxy;

Rxyest = Rxx*h;
s2 = (h'*Rxyest - 2*h'*Rxy + y'*y)'/length(y);
% s2 = mean((x*h- y).^2);
Bias = (Rxy - Rxx*h)' / Rxx * (Rxy - Rxx*h)/s2/2; 

end

