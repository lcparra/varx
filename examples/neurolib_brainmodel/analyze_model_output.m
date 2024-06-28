addpath('../../matlab/')

%%
load ../../data/neurolib_5min_rest_model_output.mat
% y = output(:,5*fs:end)'; % cut our the firt 5 seconds
y = output'; clear output

%
%y = y-mean(y);
%y = y./std(y);
model = varx(y,2);

%%
clf

h(1)=subplot(1,4,1); 
imagesc(sqrt(Cmat)); title('True connectivity')
ylabel('Brain areas'); xlabel('Brain areas');
colormap('hot')


h(2)=subplot(1,4,2);
R = model.A_Rvalue-diag(diag(model.A_Rvalue));
imagesc(R);  title('VARX estimate R'); 
xlabel('Brain areas'); 
colormap('hot')

h(3)=subplot(1,4,3);
plot(sqrt(Cmat(:)),R(:),'.'); 
xlabel('True connectivity'); ylabel('VARX estimate R')
title('Model connectivity')

h(4)=subplot(1,4,4); 
lambda = 0.01;
SC = graphicalLasso(cov(y), lambda,200); 
SC = abs(SC - diag(diag(SC))); % structural connectivity
imagesc(SC); title('Inverse covariance')
xlabel('Brain areas'); ylabel('Brain areas')
colormap('hot')

%%
set(gcf, 'CurrentAxes', h(1)); text(-10,0,'A)','FontSize',14,'FontName','Arial','FontWeight','bold');
set(gcf, 'CurrentAxes', h(2)); text(-10,0,'B)','FontSize',14,'FontName','Arial','FontWeight','bold');
set(gcf, 'CurrentAxes', h(3)); text(-0.2,0.5,'C)','FontSize',14,'FontName','Arial','FontWeight','bold');
set(gcf, 'CurrentAxes', h(4)); text(-10,0,'D)','FontSize',14,'FontName','Arial','FontWeight','bold');
%%

display(['inverse cov, Spearman r=' num2str(corr(SC(:),Cmat(:),'Type','Spearman'),2)])
display(['varx, Spearman r='        num2str(corr(R(:),Cmat(:),'Type','Spearman'),2)])

exportgraphics(gcf,'../../figures/neurolib_vs_varx.png','Resolution',300)


