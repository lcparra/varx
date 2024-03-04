function [A,B,A_pval,B_pval,A_D,B_D,T] = varx(Y,na,X,nb,lambda)

% [A,B] = varx(Y,na,X,nb,lambda) fits an vectorial ARX model to the MIMO
% system output Y with input X by minimizing the equation error e(t), i.e.
% equation error model:
%
% Y(t) = A*Y(t-1) + B*X(t) + e(t)
%
% where * represents a convolution.  Model (filter) parameters A and B are
% found with conventional least squares with ridge regression. They are
% stored as tensor of size [na,ydim,ydim] and [nb,ydim,xdim] respectively.
% na and nb are the legth of the filters. lambda is the regularization for
% the ridge (shrinkage) regularization and defaults to 0 and should not be
% selected larger than 1. Note that x(t) represents the history including
% the current sample in the input. Thus, we are allowing for instant
% effects. This is the norm in the signal processing literature but no in
% the Granger Causality VAR models, although there is no theoretical reason
% not to include instant effect in the external input. To avoid instant
% effects, the user can simply delay the input by one sample. 
% 
% [A,B,,A_pval,B_pval] = varx(Y,na ...) computes P-values for each channel
% (for all delays together) using the Deviance formalism.
%
% varx(Y,na,X,base,lambda) If base is not a scalar, it is assumed that it
% represent basis functions for filters B of size [filter length, number of
% basis functions]. B will have size [size(base,2),ydim,xdim], i.e. as many
% parameers for each path as basis functions. The actual filters can be
% obtained as tensorprod(base,B,2,1);
%
% varx(Y,na,X,nb,lambda) If Y is a cell array, then the model is fit on all
% data records in X and Y. All elements in the cell arrays X and Y have to
% have the same xdim and ydim, but may have different numer of rows (time
% samples).
%
% [A,~,A_pval] = varx(Y,na) Only fitst the AR portion. To provide lambda,
% set set x=[] and nb=0.
% 
% If the intention is to only fit a MA model, then the Granger formalism
% requires at least an AR portion for each output channel, without the
% interaction between ouput channels. If that is the intention, then one
% should call this function for each output channel separatelly, e.g.
% varx(y(:,i),na,x,nb)
% 
% [~,~,~,~,A_D,B_D,T] = varx(...) Returns the Deviance and number of sample
% used in the estimation, where D/T can serve as a measure of effect size.

% (c) July 10, 2023 Lucas C Parra
%     07/11/2023 put duplicade code into model_fit function
%     08/01/2023 transposed A, B to match varx_filter.m and varm.m, 
%     08/02/2023 corrected scaling of Bias term. 
%     08/06/2023 scalling gamma=lambda/(T-sum(lags)) for correct pvalue for large T
%     08/09/2023 error computed using Rxx, Rxy, h, and scaled LLR with T-sum(lags)+lags(i) for correct pvalues for small T and unequal lags
%     08/10/2023 moved computation of yest to varx_simulate()
%     08/16/2023 converted format of A,B to 3D tensor
%     08/18/2023 added basis functions
%     08/20/2023 removed the need to pass x,y to fit_model to save memory; added option to operate on multiple datarecords
%     08/23/2023 allow for AR portion of the model only
%     08/31/2023 changed input history to include current sample x(t) as predictor
%     09/01/2023 ignore x argument for nb=0 
%     09/11/2023 changed the order of input arguments
%     regularization range be between 0-1 with increasing T (based on
%     simulaitons)
%     12/10/2023 computation of y power now inside myxcorr() to allow removal of NaN
%     12/12/2023 determine T inside myxcorr() to account for data actually used in the fit. 
%     01/09/2024 returns D and T, with D/T as a rough measure of effect size
%     04/03/2024 Changed to Tikhonov regularization so that all variables have same regularization regardless of power. pvalues otherwise not correct when power very different after applying basis functions 

% If not simulating eXternal MA channel then xdim=0
if nargin<3 | nb==0, X=[]; nb=0; end 

% make it into a cell, if not already, so we can deal with list of data records using the same code
if ~iscell(Y), Y={Y}; end
if ~iscell(X), X={X}; end

ydim = size(Y{1},2); 
xdim = size(X{1},2); 

% generated basis tranform if using basis functions for B
if length(nb)>1   
    for i=1:ydim, base{i}=eye(na); end
    for i=1:xdim, base{end+1} = nb; end
    [nb,bparams]=size(nb); % lags according and number of parameters according to basis
else 
    base = cell(ydim+xdim,1); % empty bases
    bparams = nb; % number of parameters same as number of lags
end

% number of lags and number of parameters equal ...
lags   = ones(ydim,1)*na; params = lags;

% ... unless using basis function and need only including when modeling MA of external input
if nb, lags = [lags; ones(xdim,1)*nb]; params = [params; ones(xdim,1)*bparams]; end

tic
% calculate the correlations
Rxx=0; Rxy=0; ryy=0; T=0;
for i=1:length(Y)

    % Set preceeding output and input both as input to the LS problem
    x = Y{i}(1:end-1,:);
    y = Y{i}(2:end,:);
    
    % if modeling also the MA of eXternal input
    if nb, x = [x X{i}(2:end,:)]; end

    % Compute auto and cross correlations
    [Rxx_,Rxy_,ryy_,T_] = myxcorr(x,y,lags); 

    % accumulate over all data records
    Rxx = Rxx + Rxx_; Rxy = Rxy + Rxy_; ryy = ryy + ryy_; T = T + T_;
end
my_toc('time to compute correlations',5)


if ~exist('lambda','var') || isempty(lambda) 
    gamma = 0; % no regularization
else 
    gamma = lambda/sqrt(T-sum(lags)); % regularization decreasing with degrees of freedom
    if ~isempty(base{1})
    end
end

tic 
% fit model and produce some stats for Deviance calculation
[AB,s2,Bias] = fit_model(Rxx,Rxy,ryy,gamma,base);
my_toc('time to fit VARX model',0.5)

% store the filter values for the output as size [n,ydim,xdim]
A = permute(reshape(AB(1:ydim*na,    :),na,ydim,ydim),[1 3 2]);
B = permute(reshape(AB(1+ydim*na:end,:),params(end),xdim,ydim),[1 3 2]);

tic
% Granger Causal test for all inputs (external and recurrent)
if nargout>2 % if p-values requested
    xdim = size(x,2);
    for i=xdim:-1:1 % same as above but with reduced model removing i-th input
        ii = 1:sum(lags); ii((1:lags(i))+sum(lags(1:i-1)))=[]; % use inputs excluding i-th
        [~,s2r,Biasr] = fit_model(Rxx(ii,ii),Rxy(ii,:),ryy,gamma,{base{[1:i-1 i+1:xdim]}});
        df = T-sum(params); % degrees of freedom of the full model 
        D(:,i) = df*log(s2r./s2) - T*Biasr + T*Bias; % not the exact formula, but I calibrated and seems to work well for small T
        pval(:,i) = 1-chi2cdf(D(:,i),params(i));
    end
    A_pval = pval(:,1:ydim); 
    A_D    =    D(:,1:ydim); 
    B_pval = pval(:,ydim+1:end); 
    B_D    =    D(:,ydim+1:end); 
end
my_toc('time to compute Granger p-values',5)

function [h,s2,Bias] = fit_model(Rxx,Rxy,ryy,gamma,base)

% apply basis functions, if available 
if ~isempty(base{1})
    B = blkdiag(base{:});
    Rxx = B'*Rxx*B;
    Rxy = B'*Rxy;
end

% Regularizer
% Gamma = gamma*eye(size(Rxx))*mean(real(eig(Rxx))); % ridge regularizer, suppress extreme values
Gamma = gamma*diag(diag(Rxx)); % Tikhonov, scaled for all variables to be regularized equally, regardless of magnitude

% Least squares estimate with regularization
h = (Rxx+Gamma)\Rxy;

% mean error square
Rxyest = Rxx*h;
s2 = (sum(h.*Rxyest) - 2*sum(h.*Rxy) + ryy)';

if sum(s2<0)
    warning('Square error of linear regression can not be negative. Please debug ...')
    keyboard
end

% Bias term for ridge regression bias -- see Babadi derivation 
if gamma>0, Bias = sum((Rxy - Rxyest) .* (Rxx\(Rxy - Rxyest)))'./s2/2; else, Bias=0; end

function my_toc(message,min_time)

t=toc; if t>min_time, disp([message ': ' num2str(t)]); end



