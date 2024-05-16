function [Rxx,Rxy,ryy,T] = myxcorr(x,y,lags)
% [Rxx,Rxy,ryy,T] = myxcorr(x,y,lags) 
% 
% Computes auto- and cross-correlation matrices Rxx, Rxy and ryy defined as
% (after the mean is subtracted from x and y):
%
%   Rxx(l)= sum_n x'(n)*x(n+l)
%   Rxy(l)= sum_n x'(n)*y(n+l)
%   ryy   = sum_n |y(n)|^2
%
% This is computed in the time domain, returning only zero and
% positive delays l = 0 ... lags-1. Note that these definitions differ from
% the standard definition of the cross correlation and from matlab's
% conventional xcorr() implementation. Here y is delayed, not x, which is
% equivalent to keeping the negative delays in the standard definition.
%
% Rxx and Rxy are arranged as required to implement conventional MIMO
% system identification (see example below). Rxx is a square block-toeplitz
% auto-correlation matrix of size sum(lags) and Rxy is a cross-correlation
% block-matrix of size [sum(lags), ydim]. lags is a vector of lags used for
% each of xdim columns in x. 
%
% x and y are required input, all others are optional. lags defaults to
% size(x,1)-1 and is converted into a vector of size xdim if it is given as
% a scalar.  x and y should be arranged as time by dimensions. y should
% not be longer in time than x; if necessary, pad x with zeros before
% calling.
%
% Correlations are computed using Toeplitz marices. x or y may contain NaN
% which are removed from the sum over n which have NaN in any row of y, or
% do not have a max(lags) history in the rows of x. T is returned to report
% how many samples were used. The max(lags) samples at the start are also
% removed (this is known as the covariance method, see Ljung, System
% Identification, Ch. 10.1). Mean subtraction is done with the mean of
% valid samples only.
%
% This functon can be used for MIMO FIR identification as follows: 
% 
% [Rxx,Rxy] = myxcorr(x,y,lags);
% h = inv(Rxx)*Ry;
% yest = filterMIMO(h,x,lags);
% error = y-yest;
 
% (c) April 21, 2021, Lucas C Parra
%     04/27/2021, added variable lags for each x dimension
%     05/04/2021, squeezed to improove speed and memory 
%     06/14/2021, fixed a bug with symmetry of Rxx
%     07/09/2023, moved expand function outside, so I can generalize for use in Granger causality estimation
%     08/30/2023, allow for scalar lag argument
%     12/10/2023, compute sum of squares of y, and handle NaN with Toeplitz
%     12/11/2023, subtract the mean of x,y prior to computing anything. 
%     12/12/2023, removed fft implementation, only using Toeplitz, omit samples at start  
  

if ~exist('lags','var') || isempty(lags), lags=size(x,1)-1; 
elseif length(lags)==1, lags=lags*ones(size(x,2),1); end

% will compute correlations up to the largest possible lag
Q=max(lags);

% valid samples means that there are no NaN in all rows of Y and Q history of all X rows
z = nan(Q-1,1); []; % removed starting values (with is what varm() does in matlab. 
valid = ~isnan(filter(ones(Q,1),1,sum(x,2),z)+sum(y,2)); 

% Number of valid samples
T = sum(valid);

% remove offset
x = x - mean(x(valid,:));
y = y - mean(y(valid,:));

% compute correlations with block-toplitz matrix
for i=size(x,2):-1:1
    X(:,(1:lags(i))+sum(lags(1:i-1))) = toeplitz(x(:,i),[x(1,i) zeros(1,lags(i)-1)]);
end
Rxx = X(valid,:)'*X(valid,:);
Rxy = X(valid,:)'*y(valid,:);

% Power of y
ryy = sum(abs(y(valid,:)).^2);

return

%% debuging code
h = Rxx\Rxy;
Rxyest = Rxx*h;
s2 = (sum(h.*Rxyest) - 2*sum(h.*Rxy) + ryy)';
if sum(s2<0)
    keyboard
end


%% test code for MIMO MA identification 
clear all
load SteamEng % arbitrary data to demo the code, not to test performance!
x=[Pressure,MagVolt]; % some multidimentional input
y=[GenVolt,Speed]; % x(30,1)=NaN;

% Least squares MIMO - MA estimate with regularization
L = 10; % FIR filter order
[Rxx,Rxy,ryy,T] = myxcorr(x,y,L); % compute correlations
D = eye(size(Rxx))*mean(eig(Rxx))*0; % 0.5 shrinkage regularization
h = (Rxx + D)\Rxy; % regularized LS estimate
[yest,H] = filterMIMO(h,x,L);
SS = sum(nansum((yest-y).^2))

% show some results
clf
ydim = size(y,2);
for i=1:ydim
    subplot(3,size(y,2),1); 
    plot(x) % compare original to estimated output
    title('Inputs');    xlabel('time (s)')
    subplot(3,ydim,i+ydim*1); 
    plot([y(:,i) yest(:,i)]) % compare original to estimated output
    title(['Target and estimated output ' num2str(i)])
    subplot(3,ydim,i+ydim*2); 
    stem([H{i,:}]) % show estimated filters
    title(['FIR filter estimates ' num2str(i)])
    xlabel('time delay (s)')
end

