clear all
addpath('../matlab/')

% Here we define a 4 node network with 1&2 intereacting separately from 3&4
A(:,1:2,1) = [[0.3 -0.5]',[0.1 0 ]']; 
A(:,1:2,2) = [[-.5 0.4]',[0.5 -.7 ]'];
A(:,3:4,3) = [[0.3 -0.5 ]',[0.1 0 ]']; 
A(:,3:4,4) = [[-.5 0.4 ]',[0.5 -.7 ]'];
% a single input only affects node 2
B = [[0 0 0]',[1 0.5 0.2]',[0 0 0]',[0 0 0]']; 

% some dimensions 
[nb,ydim,xdim] = size(B);
[na,ydim,ydim] = size(A);
T = 10000;
x = randn(T,xdim);

% implement gain adaptation in time window of 1/gamma samples 
gamma = 0.001; 

% simulate with stimulus and estimate the model
y_stim = varx_simulate(B,A,x,1,gamma);  
model_stim = varx(y_stim,na,x,nb);
varx_display(model_stim)

% simulate "rest" without stimulus and estimate the model
y_rest = varx_simulate(B,A,zeros(size(x)),1,gamma); 
model_rest = varx(y_rest,na,x,nb);

% Compute change in noise level
delta = db(model_stim.s2'./var(y_stim)) - db(model_rest.s2'./var(y_stim));

%%
nexttile(4)
bar(delta); hold on
ax=axis; plot(ax(1:2),[0 0],'k')
title('e^2/y^2 power ratio change')
xlabel('Channels'); ylabel('dB'); 
delete(nexttile(7)) 

% grab the handles of all axes 
h=get(get(gcf,'Children'),'Children');
set(gcf, 'CurrentAxes', h(5)); text(-2,0,'A)','FontSize',14,'FontName','Arial','FontWeight','bold');
set(gcf, 'CurrentAxes', h(3)); text(-3,0,'B)','FontSize',14,'FontName','Arial','FontWeight','bold');
set(gcf, 'CurrentAxes', h(1)); text(-0.7,1,'C)','FontSize',14,'FontName','Arial','FontWeight','bold');

exportgraphics(gcf,'../figures/simulate_gain_adaptation.png','Resolution',300)


% Finding: the noise estimates are unchanged, but the power is higher with
% external input, so the innovation noise power is relatively less with
% external input when using gain adaptation to maintain the total power
% constant. If the gain adaptation is done in real time, then the noise
% reduction stays localize with the channel that has external input. If the
% adaptation is not after the fact, then the power of the input spreads
% through the network, and all nodes see a reduction of noise. Interesting
% prediction for the neural data. 
