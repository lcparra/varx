addpath('../../matlab/')

%%
load ../../data/neurolib_5min_rest_model_output.mat
y = output(:,5*fs:end)'; % cut our the firt 5 seconds

%
y = y-mean(y);
y = y./std(y);
model = varx(y,2);
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
plot(R(:),sqrt(Cmat(:)),'.'); 
xlabel('True connectivity'); ylabel('VARX estimate R')
title('Model connectivity')

h(4)=subplot(1,4,4); 
lambda = 0.01;
SC = graphicalLasso(cov(y), lambda,200); 
SC = abs(SC - diag(diag(SC))); % structural connectivity
imagesc(SC); title('Inverse covariance')
xlabel('Brain areas'); ylabel('Brain areas')
colormap('hot')


sublabel(h,-5,-25);

display(['inverse cov, Spearman r=' num2str(corr(SC(:),Cmat(:),'Type','Spearman'),2)])
display(['varx, Spearman r='        num2str(corr(R(:),Cmat(:),'Type','Spearman'),2)])

saveas(gca,'../../figures/neurolib_vs_varx.png')



