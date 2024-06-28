function pval = granger_test(x,y,Q,z,control)
% pval = granger_test(x,y,Q,z,control)
% 
% Test for significance of a Granger "causal" link, from X to Y. In other
% words, it tests if time series Y can be linearly predicted from Q
% preceeding values in time series X, above and beyond what Q past values
% of Y predict. The test statistics is the Deviance, i.e. the
% log-likelihood ratio which follows as Chi-square distribution, with which
% the pvalue is computed. One can additionally "control" for a covariate z,
% in one of three ways, specified in control: 'instant' adds the current
% and Q-1 preceeding values of z as predictor. 'delayed' adds the Q
% preceeding values of z as predictor. 'none' does nothing and ignores z.

% (c) Lucas Parra, May 9, 2022. All errors mine. Use at your own risk.

if ~exist('control'), control='none'; end

X = toeplitz(x(Q:end-1),x(Q:-1:1));
Y = toeplitz(y(Q:end-1),y(Q:-1:1));
switch control
    case 'instant'
        Z = toeplitz(z(Q+1:end),z(Q+1:-1:2));
    case 'delayed'
        Z = toeplitz(z(Q:end-1),z(Q:-1:1));
    case 'none'
        Z = []; %
    otherwise
        error('Control must be set to one of these: instant, delayed, none.')
end
bf = [X Y Z]\y(Q+1:end); % full model
br = [  Y Z]\y(Q+1:end); % reduced model
SSf = sum((y(Q+1:end)-[X Y Z]*bf).^2);
SSr = sum((y(Q+1:end)-[  Y Z]*br).^2);
[T,P] = size(X); % making sure numbers are correct for stats
D = T*log(SSr/SSf); % Deviance (2*log-likelihood ratio)
pval = 1-chi2cdf(D,P);

end