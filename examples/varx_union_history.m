addpath('../matlab/')
clear all
warning off % warnings about headers are annoying me
Tab = readtable('../data/union-history-US.xlsx');
warning on
data = Tab.Variables;

% select variables
year = data(:,1);
Y = data(:,[2 3 4]);
ynames = {Tab.Properties.VariableNames{[2 3 4]}};
X = data(:,5);
xnames = {Tab.Properties.VariableNames{5}};

% One way of choosing order of the model; there are others
%[pacf,lags,bounds] = parcorr(Y(:,3),10);
%na = max(lags((bounds(1)<pacf | pacf<bounds(2))));

na=3;
nb=na; % results too dependent on order, need larger dataset. 
lambda=0.00; % I don't like that causality direction depends on this parameter
model = varx(Y,na,X,nb,lambda);

figure(1)
[Graph,G_plot]=varx_display(model,yname=ynames,xname=xnames,plottype='Graph',threshold=0.05,duration=10);

figure(2)
clf
h(1)=subplot(1,2,1); plot(year,[Y X]); stkplt
legend({ynames{:},xnames{:}},'Location','southoutside')
xlabel('Year')
axis tight
h(2)=subplot(1,2,2);
plot(Graph,'LineWidth',G_plot.width*0.5,'XData',G_plot.xdata,'YData',G_plot.ydata,...
            'EdgeColor',G_plot.color,'NodeColor',G_plot.nodecolor);
axis off

sublabel(h,10,-20);
exportgraphics(gcf,'../figures/varx_union_history.png','Resolution',300)