addpath('../../matlab/')

%%
load ../../data/neurolib_model_output_rest.mat
% load ../../data/neurolib_model_output_asymmetry_0.10_rest.mat
% y = output(:,5*fs:end)'; % cut our the firt 5 seconds
y = output'; clear output

%
%y = y-mean(y);
%y = y./std(y);
model = varx(y,2);

%% graphical lasso to do sparce inverse covariance:
lambda = 0.01;
invC = graphicalLasso(cov(y), lambda,200); 
invC = abs(invC - diag(diag(invC))); % structural connectivity


%%
figure(1)
clf

load('colormap_reds.mat','colormap_reds')

h(1)=subplot(2,4,1); 
imagesc(sqrt(Cmat)); title('Struct. Connectivty C')
clim([0 1])
ylabel('Brain areas'); xlabel('Brain areas');
colormap(gca,colormap_reds)


h(2)=subplot(2,4,2);
R = model.A_Rvalue-diag(diag(model.A_Rvalue));
imagesc(R);  title('VARX estimate R'); 
clim([0 0.4])
xlabel('Brain areas'); 
colormap(gca,colormap_reds)

h(3)=subplot(2,4,3);
plot(sqrt(Cmat(:)),R(:),'.'); 
xlim([0 1]); ylim([0 0.4])
xlabel('True C'); ylabel('Model R')
title('Validation')

h(4)=subplot(2,4,4); 
imagesc(invC); title('Sparce Inverse FC')
xlabel('Brain areas'); ylabel('Brain areas')
colormap(gca,colormap_reds)

h(5)=subplot(2,4,5); 
imagesc(sqrt(Cmat)-sqrt(Cmat)'); title('C - C^T')
clim([-1 1]*0.8)
ylabel('Brain areas'); xlabel('Brain areas');
colormap(gca,redblue)


h(6)=subplot(2,4,6);
R = model.A_Rvalue-diag(diag(model.A_Rvalue));
imagesc(R-R');  title('R-R^T'); 
clim([-1 1]*0.3)
xlabel('Brain areas'); 
colormap(gca,redblue)

h(7)=subplot(2,4,7);
Cmat_diff=Cmat-Cmat'; R_diff=R-R';
plot(Cmat_diff(:),R_diff(:),'.'); 
xlabel('C-C^T'); ylabel('R - R^T')
title('Validation')

h(8)=subplot(2,4,8); 
imagesc(invC-invC'); title('inv FC - inv FC^T')
xlabel('Brain areas'); ylabel('Brain areas')
colormap(gca,redblue)

%
set(gcf, 'CurrentAxes', h(1)); text(-20,0,'A)','FontSize',14,'FontName','Arial','FontWeight','bold');
set(gcf, 'CurrentAxes', h(2)); text(-20,0,'B)','FontSize',14,'FontName','Arial','FontWeight','bold');
set(gcf, 'CurrentAxes', h(3)); ax=axis; text(ax(1)-0.25,ax(4),'C)','FontSize',14,'FontName','Arial','FontWeight','bold');
set(gcf, 'CurrentAxes', h(4)); text(-20,0,'D)','FontSize',14,'FontName','Arial','FontWeight','bold');

set(gcf, 'CurrentAxes', h(5)); text(-20,0,'E)','FontSize',14,'FontName','Arial','FontWeight','bold');
set(gcf, 'CurrentAxes', h(6)); text(-20,0,'F)','FontSize',14,'FontName','Arial','FontWeight','bold');
set(gcf, 'CurrentAxes', h(7)); ax=axis; text(ax(1)-0.35,ax(4),'G)','FontSize',14,'FontName','Arial','FontWeight','bold');
set(gcf, 'CurrentAxes', h(8)); text(-20,0,'H)','FontSize',14,'FontName','Arial','FontWeight','bold');

%
display(['inverse cov, Spearman r=' num2str(corr(invC(:),Cmat(:),'Type','Spearman'),2)])
display(['varx, Spearman r='        num2str(corr(R(:),Cmat(:),'Type','Spearman'),2)])

% exportgraphics(gcf,'../../figures/neurolib_vs_varx.png','Resolution',300)
% exportgraphics(gcf,'../../figures/neurolib_vs_varx_asymmetry.png','Resolution',300)


%% Show the A matrix

figure(2)
h(1)=subplot(1,3,1);
R = model.A_Rvalue;
imagesc(R);  title('VARX estimate R'); 
clim([0 0.4])
ylabel('Brain areas'); 
xlabel('Brain areas'); 
colormap(gca,colormap_reds)
colorbar

h(2)=subplot(1,3,2);
imagesc(squeeze(model.A(1,:,:)));  title('VARX A(1)'); 
xlabel('Brain areas'); 
clim([-1 1]/4)
colormap(gca,'parula')
colorbar

h(3)=subplot(1,3,3);
R = model.A_Rvalue-diag(diag(model.A_Rvalue));
imagesc(squeeze(model.A(2,:,:)));  title('VARX A(2)');  
xlabel('Brain areas'); 
clim([-1 1]/4)
colormap(gca,'parula')
colorbar
exportgraphics(gcf,'../../figures/neurolib_A_matrix.png','Resolution',300)
