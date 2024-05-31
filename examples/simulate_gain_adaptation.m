clear all

% Here we define a 4 node network with 1&2 intereacting separately from 3&4
A(:,1:2,1) = [[0.3 -0.5 0]',[0.1 0 0]']; 
A(:,1:2,2) = [[-.5 0.4 0]',[0.5 -.7 0]'];
A(:,3:4,3) = [[0.3 -0.5 0]',[0.1 0 0]']; 
A(:,3:4,4) = [[-.5 0.4 0]',[0.5 -.7 0]'];
% a single input only affects node 2
B = [[0 0]',[1 1]',[0 0]',[0 0]']; 

% some dimensions 
[nb,ydim,xdim] = size(B);
[na,ydim,ydim] = size(A);
T = 10000;

% implement gain adaptation in time window of 1/gamma samples 
gamma = 0.001; 

% simulate with stimulus and estimate the model
x = randn(T,xdim);
y = varx_simulate(B,A,x,1,gamma);  
model_stim = varx(y,na,x,nb);

% simulate "rest" without stimulus and estimate the model
y = varx_simulate(B,A,zeros(T,xdim),1,gamma); 
model_rest = varx(y,na,x,nb);

% Compute change in noise level
model_stim.s2 - model_rest.s2

% Finding: the noise estimates are unchanged, but the power is higher with
% external input, so the innovation noise power is relatively less with
% external input when using gain adaptation to maintain the total power
% constant. If the gain adaptation is done in real time, then the noise
% reduction stays localize with the channel that has external input. If the
% adaptation is not after the fact, then the power of the input spreads
% through the network, and all nodes see a reduction of noise. Interesting
% prediction for the neural data. 
