clf
load Data_USEconModel
valid = ~isnan(sum(Data,2));
Time =  DataTimeTable.Time(valid);
sort = [1 2 4 5 6 7 9 11 12 3 8 14]; % rates last, not M1 money supply 
Data= Data(valid,sort);
series = {series{sort}};
N = size(Data,2);

% get the abraviations
for i=1:length(series), 
    s=sscanf(series{i},'(%s'); 
    short_name{i}=s(1:end-1); 
end; 

% convert all into anual rates
%Data = [Data(5:end,1:10)./Data(1:end-4,1:10) Data(5:end,11:14)];
%Time = Time(5:end);

% remove 2008 and 2009
%Data = Data(1:end-8,:);
%Time = Time(1:end-8);


addpath('../matlab/')
exogenous = [3 10]; % Goverment Expenditure; Federal Funds Rate
endogenous = setdiff(1:N,exogenous); 


Y = Data(:,endogenous); yname={short_name{endogenous}};
X = Data(:,exogenous ); xname={short_name{exogenous}};
model = varx(Y,6,X,6);

% some displays
figure(1)
[Graph,G_plot]=varx_display(model,plottype='Graph no filter',xname=xname,yname=yname,threshold=0.001)

clf
h(1)=subplot(2,1,1); %tiledlayout(1,3); nexttile([1 1])
semilogy(Time,Data); xlabel('year')
ylim([10^-1 10^4])
legend(series','Location','eastoutside','box','off')
title('History')

% nexttile([1 2])
h(2)=subplot(2,2,3); 
plot(Graph,'LineWidth',G_plot.width,'XData',G_plot.xdata,'YData',G_plot.ydata,...
             'EdgeColor',G_plot.color,'NodeColor',G_plot.nodecolor);
ylim([-1 1]); axis off; title('Effects')

h(3)=subplot(2,2,4); 
imagesc(corrcoef(Data)); colorbar; axis square
title('Correlation'); set(gca,'XTickLabel',short_name); set(gca,'YTickLabel',short_name)
set(gca,'XTick',1:14);  set(gca,'YTick',1:14); 
set(gca,'XTickLabelRotation',45); set(gca,'FontSize',8)


sublabel(h(1:3),-10,-40)


exportgraphics(gcf,'../figures/economy_example.png', 'Resolution', 300)


% for i=1:N, text(2,1.5-i*0.2,series{i}); end; 