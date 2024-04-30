
%% VARX identification on real data
clear all
load SteamEng % arbitrary data to demo the code, not to test performance!
x=[Pressure,MagVolt]; xname = {'Pressure','MagVolt'}; 
y=[GenVolt,Speed]; yname = {'GenVolt','Speed'};
na=10; nb=20; lambda=0;

model = varx(y,na,x,nb,lambda);

yest = varx_simulate(model.B,model.A,x,y);

figure(1); show_prediction(x,y,yest);
figure(2); varx_display(model,xname=xname,yname=yname,duration=nb,plottype='Graph');


%% generate a data with known VARX model and compare varx result with matlab's varm estimate
clear all

% define VARX model. Where there are all zeros, that path does not
% contribute
A(:,:,1) = [[0.3 -0.5 0]',[0 0 0]']; A(:,:,2) = [[-.5 0.4 0]',[0.5 -.7 0]'];
B = [[1 0]',[0 0]']; % single input. the varm estimate only works for zero delay
[nb,ydim,xdim] = size(B);
[na,ydim,ydim] = size(A);

% simulate 
lambda = 0.0;
T = 1000;
x = randn(T,xdim);
[y,e] = varx_simulate(B,A,x,1);  

% estimate VARX model
model = varx(y,na,x,nb,lambda);

[yest,eest] = varx_simulate(model.B,model.A,x,y);

R2=1-sum(eest.^2)./sum(y.^2)

% estimate model using Econometric toolbox' mvar
ydim = size(A,2);
Mdl = varm(ydim,na);
EstMdl = estimate(Mdl,y,'X',x);

% some displays
figure(3); varx_display(model,plottype='matrix',xname={'x1'},yname={'y1','y2'});
exportgraphics(gcf,'../figures/known-model-efficacy.png', 'Resolution', 300)
figure(4); show_prediction(x,y,yest);

% compare estimate to truth
figure(2); clf
for t=1:na, for i=1:ydim, for j=1:ydim, AR{i,j}(t) = EstMdl.AR{t}(i,j); end; end; end;
plot(A(:),model.A(:),'.'); hold on;
plot(A(:),[AR{:}],'.');
plot(B(:),model.B(:),'x')
plot(B(1,:),EstMdl.Beta(:),'x')
plot([-1 1],[-1 1]); hold off
axis equal; axis tight
legend('our {\bf A}','matlab {\bf A}','our {\bf B}','matlab {\bf B}','Location','eastoutside')
xlabel('true value'); ylabel('estimate')
exportgraphics(gca,'../figures/known-model-true-vs-estimate.png', 'Resolution', 300)

%% determine the effect of regularization
clear A B
A(:,:,1) = [[0.8 0.1]',[-0.5 0.2]']; A(:,:,2) = [[0 0]',[-0.3 0.2]'];
B(:,:,1) = [ones(2,1),[1;zeros(1,1)]]; % single input. the varm estimate only works for zero delay
B(:,:,2) = [ones(2,1),zeros(2,1)]; % single input. the varm estimate only works for zero delay
[nb,ydim,xdim] = size(B);
[na,ydim,ydim] = size(A);

T = [30 50 100];
lambda = [0:0.2:1];

% simulate 
E_test=zeros(length(lambda),length(T),1000);
E_train=zeros(length(lambda),length(T),1000);
clear Amissing Apresent Bmissing Bpresent
for nt = 1:length(T)
    for nrand = 1:size(E_train,3)
        x = randn(2*T(nt),xdim);
        [y,e] = varx_simulate(B,A,x,2); 
  %      x(30,1) = NaN;         % test if addeind nan works 
        for i=1:length(lambda)
            model = varx(y(1:T(nt),:),na,x(1:T(nt),:),nb,lambda(i));
            Apresent(i,nt,nrand)=model.A_pval(1,2); 
            Amissing(i,nt,nrand)=model.A_pval(2,1);
            Bpresent(i,nt,nrand)=model.B_pval(2,2); 
            Bmissing(i,nt,nrand)=model.B_pval(1,2);
            [yest,eest] = varx_simulate(model.B,model.A,x,y);
            E_train(i,nt,nrand)= mean(sum(eest(1:T(nt)    ,:).^2)./sum(y(1:T(nt)    ,:).^2));
            E_test (i,nt,nrand)= mean(sum(eest(T(nt)+1:2*T(nt),:).^2)./sum(y(T(nt)+1:2*T(nt),:).^2));
        end
    end
    clear str; for i=1:length(T), str{i}=['T=' num2str(T(i))]; end;
    subplot(2,2,1)
    errorbar(lambda, mean(E_train,3),std(E_train,[],3)/sqrt(size(E_train,3))); 
    ylabel('Training Error'); xlabel('\lambda')
    legend(str{:},'location','northwest')
    subplot(2,2,2)
    errorbar(lambda, mean(E_test, 3),std(E_test ,[],3)/sqrt(size(E_test ,3))); 
    ylabel('Test Error'); 
    xlabel('\lambda')
    subplot(2,2,3); 
    plot(lambda,mean(Apresent<0.05,3)); hold on
    plot(lambda,mean(Bpresent<0.05,3),':'); 
    plot([min(lambda) max(lambda)],[1 1]*0.05,'--k'); hold off
    ylim([0 0.2]); xlabel('\lambda'); ylabel('False discovery rate')
    subplot(2,2,4); 
    plot(lambda,mean(Amissing<0.05,3)); hold on 
    plot(lambda,mean(Bmissing<0.05,3),':'); hold off
    ylim([0 1.1]); xlabel('\lambda'); ylabel('Statistical power')
    drawnow
end

exportgraphics(gcf,'../figures/effect_of_regularization.png','Resolution',300)


%% try repeating this 1000 times and see if the pvalues are correct
clear all
A(:,:,1) = [[0.9 -0.5 0]',[0 0 0]']; A(:,:,2) = [[-.5 0.4 0]',[0.5 -.7 0]'];
B = [[1 0 0 0 -1 0 0 0 0 1 0 0 -1 ]',[0 0 0 0 0 0 0 0 0 0 0 0 0 ]']; % the varm estimate only works for zero delay
[nb,ydim,xdim] = size(B);
[na,ydim,ydim] = size(A);
clear A_pval B_pval
lambda = 1;
for n=1000:-1:1
    x = randn(100,size(B,3)); 
    y = varx_simulate(B,A,x,1); 
    x(30,1) = NaN; % test if adding NaN works
    model = varx(y,na,x,nb,lambda);
    A_pval(:,:,n)=model.A_pval; B_pval(:,:,n)=model.B_pval;
end
mean(A_pval<0.05,3)
mean(B_pval<0.05,3)

%% now do it with some larger models
tiledlayout(1,6);
na = 2; nb = 2; xdim=1;
for ydim=[6 60];
    clear A B
    for i=1:ydim,
        for j=1:ydim, A(:,i,j) = 0.05*sign(randn(na,1)); end;
        for j=1:xdim, B(:,i,j) = randn(nb,1); end;
    end
    % set two of the paths to zero, to see if we correctly identify them at alpha<0.05
    A(:,2,2)=0;
    B(:,5,1)=0;
    lambda=0;
    clear A_pval B_pval
    for n=1000:-1:1
        x = randn(1000,xdim);
        y = varx_simulate(B,A,x,1);
        model = varx(y,na,x,nb,lambda);
        A_pval(:,:,n)=model.A_pval; B_pval(:,:,n)=model.B_pval;
    end
    nexttile([1 2]); imagesc(mean(A_pval<0.05,3)); axis equal; axis tight; clim([0 0.1]);
    title('A')
    h = colorbar; ylabel(h,'false discovery rate at p<0.05')
    nexttile([1 1]); imagesc(mean(B_pval<0.05,3)); axis equal; axis tight; clim([0 0.1]);
    title('B'); 
    drawnow
end

exportgraphics(gcf,'../figures/false_alarm_ydim-6-60.png','Resolution',300)

%% now try using basis functions
clear A B
A(:,:,1) = [[0.9 -0.5 0]',[0 0 0]']; A(:,:,2) = [[-.5 0.4 0]',[0.5 -.7 0]'];
B = [[1 0.5 0 -0.5]',[0 0 0 0]']; % basis coefficients
[nb,ydim,xdim] = size(B);
[na,ydim,ydim] = size(A);
base = basis(30,nb,'normal');
lambda = 0.1;

x = randn(100,xdim);
y = varx_simulate(tensorprod(base,B,2,1),A,x,1);
model = varx(y,na,x,base,lambda);

yest = varx_simulate(model.B,model.A,x,y);
% x(30,1) = NaN; % test if adding NaN works


figure(3); varx_display(model);
figure(4); show_prediction(x,y,yest);

figure(5); clf
plot(A(:),model.A(:),'o'); hold on
plot(B(:),model.B_coeff(:),'o'); 
plot([-1 1],[-1 1],':k'); hold off
xlabel('True values'); ylabel('Estimated values')
legend('AR coeeficients','MA basis coefficients','location','northwest')

%% try repeating this 1000 times and see if the pvalues are correct
clear all
A(:,:,1) = [[0.9 -0.5 0]',[0 0 0]']; A(:,:,2) = [[-.5 0.4 0]',[0.5 -.7 0]'];
B = [[1 0.5 0 -0.5]',[0 0 0 0]']; % basis coefficients
[nb,ydim,xdim] = size(B);
[na,ydim,ydim] = size(A);
base = basis(60,nb,'normal');
B = tensorprod(base,B,2,1);

for n=5000:-1:1
    x = randn(300,xdim);
    % generate using basis function
    y = varx_simulate(B,A,x,1);
    % estimate using basis functions
    model_b = varx(y,na,x,base); 
    A_pval_b(:,:,n)=model_b.A_pval; B_pval_b(:,:,n)=model_b.B_pval;
    % estimate all parameters
    model = varx(y,na,x,60);
    A_pval(:,:,n)=model.A_pval; B_pval(:,:,n)=model.B_pval;
end

%%
clf
haxis(1)=subplot(1,3,1); plot(base); 
xlabel('lag (samples)'); ylabel('basis functions')

haxis(2)=subplot(1,3,2);
for i=1:2
    plot(model.B(:,i),'b'); hold on;
    plot(model_b.B(:,i),'r');
    plot(B(:,i),'g'); 
end
legend('free estimate','basis estimate','true filter')
hold off
xlabel('lag (samples)')
ylabel('filter B')

haxis(3)=subplot(1,3,3);
FDrate(1,1) = mean(A_pval(2,1,:)<0.05,3);
FDrate(1,2) = mean(B_pval(2,1,:)<0.05,3);
FDrate(2,1) = mean(A_pval_b(2,1,:)<0.05,3);
FDrate(2,2) = mean(B_pval_b(2,1,:)<0.05,3);
bar(FDrate); hold on
ylabel('False discovery rate'); ylim([0 0.07])
ax=axis; plot(ax(1:2),[0.05 0.05],':k'); hold off
legend('A link','B link','\alpha=0.05','location','northwest')
set(gca,'XTickLabel',{'free','basis'})
xlabel('Estimate')

sublabel(haxis,-5,-35);
saveas(gcf,'../figures/effect_of_basis.png')

%% demo modeling just the AR part
clear all
A(:,:,1) = [[0 0 0]',[0. -0.5 0]']; A(:,:,2) = [[-.5 0.4 0]',[0.5 -.7 0]'];
ydim = size(A,2); xdim=0;
B=zeros(0,ydim,0);
na=size(A,1); nb=size(B,1); lambda=0.5;
x=randn(200,xdim); % no external input
[y,e] = varx_simulate(B,A,x,1);

model = varx(y,na);
model.A_pval

figure(3); varx_display(model,yname={'y1','y2'},duration=20,plottype='Graph');


figure(5); clf
plot(A(:),model.A(:),'o'); hold on
plot(B(:),model.B(:),'o'); 
plot([-1 1],[-1 1],':k'); hold off
xlabel('True values'); ylabel('Estimated values')
legend('AR coeeficients','location','northwest')


%% demo modeling just the MA part
clear all
% Notice that Granger formalism still includes an AR portion for each
% channel. The following A will color the output but without coupling the output channels.
A(:,:,1) = [[0.9 -0.5 0.2]',[0 0 0]']; A(:,:,2) = [[0 0 0]',[0.5 -.7 0.1]'];
B(:,:,1) = [[1 0.5 0 -0.5]',[0 0 0 0]']; B(:,:,2) = [[0.7 0.2 0.1 -0.4]',[1 0 0 0]']; 
[nb,ydim,xdim] = size(B);
[na,ydim,ydim] = size(A);
x=randn(200,xdim); % no external input
[y,e] = varx_simulate(B,A,x,1);


model.A_pval = nan(ydim,ydim); % off dianoal elements are not computed, make sure that is clear. 
model.A_Deviance = nan(ydim,ydim); % off dianoal elements are not computed, make sure that is clear. 
for i=1:ydim
    m = varx(y(:,i),na,x,nb);
    model.A(:,i,i)=m.A;
    model.B(:,i,:)=m.B;
    model.A_pval(i,i)=m.A_pval;
    model.B_pval(i,:)=m.B_pval;
    model.A_Deviance(i,i)=m.A_Deviance;
    model.B_Deviance(i,:)=m.B_Deviance;
    model.T = m.T;
end
model.A_pval,model.B_pval

[yest,e] = varx_simulate(model.B,model.A,x,y);

figure(3); varx_display(model,yname={'y1','y2'},xname={'x1','x2'},duration=20,plottype='Graph');
figure(4); show_prediction(x,y,yest);

figure(5); clf
plot(A(:),model.A(:),'o'); hold on
plot(B(:),model.B(:),'o'); 
plot([-1 1],[-1 1],':k'); hold off
xlabel('True values'); ylabel('Estimated values')
legend('MA coeeficients','location','northwest')


%% now try fitting model with a missing input variable 

% define VARX model with no cross-talk between the two variables
clear A
A(:,:,1) = [[0.9 -0.5 0]',[0 0 0]']; A(:,:,2) = [[0 0 0]',[0.5 -.7 0]'];
B = [[1 .1]',[1 0.1]']; % single input, same effect on both output variables
%B = [[0 0]',[0 0]']; % single input, same effect on both output variables
[nb,ydim,xdim] = size(B);
[na,ydim,ydim] = size(A);

% simulate 
x = randn(1000,xdim);
[y,e] = varx_simulate(B,A,x,1);

% estimate VARX model misspecified assuming there is no input
model = varx(y,na); 

% simulate varx without the external imput (I know, its cumbersome) 
yest = varx_simulate(zeros(0,ydim,0),model.A,zeros(size(y,1),0),y);

% some displays
figure(1)
varx_display(model);
sgtitle('Simulated VARX and estimation with missing input')

% compare estimate to truth
figure(2); clf
plot(A(:),model.A(:),'.'); hold on;
plot([-1 1],[-1 1]); hold off
xlabel('true value'); ylabel('estimate')

% estimate VARX model now with the input included
model = varx(y,na,x,nb); 

% some displays
figure(3)
varx_display(model);
sgtitle('Simulated VARX and estimation with input included')

% compare estimate to truth
figure(4); clf
plot(A(:),model.A(:),'.'); hold on;
plot([-1 1],[-1 1]); hold off
xlabel('true value'); ylabel('estimate')


%% model long smooth AR coefficients, see if we recover them correctly. 
% make smooth A
clear A
A(:,:,1) = [[0.1 -0.5 0]',[.1 0 0]']; A(:,:,2) = [[-0.1 0 0]',[0.05 -.7 0]'];
[na,ydim,ydim] = size(A);
base = basis(30,na,'normal');
A = tensorprod(base,A,2,1);
[na,ydim,ydim] = size(A);

% make dummy input B, x
xdim=0;
B=zeros(0,ydim,0);
na=size(A,1); nb=size(B,1);
x=zeros(2000,xdim); % no external input

% simulate 
[y,e] = varx_simulate(B,A,x,1);

% estimate VARX model 
lambda = 0.0;
model = varx(y,na,[],[],lambda); 

% simulate equaiton error model with the estimated Aest
yest = varx_simulate(B,model.A,x,y);

% some displays
figure(1); varx_display(model);
figure(2); show_prediction(x,y,yest);

% compare estimate to truth
figure(3); clf
plot(A(:),model.A(:),'.'); hold on;
plot([-1 1],[-1 1]); hold off
xlabel('true value'); ylabel('estimate')
sgtitle('Simulated VARX and estimation')

%% Test what happens if innovation changes in power across repeats
A(:,:,1) = [[0.3 -0.5 0]',[0 0 0]']; A(:,:,2) = [[-.5 0.4 0]',[0.5 -.7 0]'];
B = [[0 0]',[1 0]']; % single input. the varm estimate only works for zero delay
[nb,ydim,xdim] = size(B);
[na,ydim,ydim] = size(A);

% simulate 
lambda = 0.0;
T = 100000;
clear x y
x{1} = randn(T,xdim); y{1} = varx_simulate(B,A,x{1}*10,1);  
x{2} = randn(T,xdim); y{2} = varx_simulate(B,A,x{2},10);  

% estimate VARX model
model = varx(y,na,x,nb,lambda);

varx_display(model,plottype='graph',xname={'x1'},yname={'y1','y2'});


% ----------------- result display function --------------------------
function show_prediction(x,y,yest)

clf
ydim = size(y,2);
xdim = size(x,2);
for i=1:ydim
    subplot(2,2,1); plot(x); title('External Inputs');
    subplot(2,2,2); plot(y); title('Recursive Input');
    subplot(2,ydim,ydim*1+i);
    plot([y(:,i) yest(:,i)]) % compare original to estimated output
    title(['Target and estimated output ' num2str(i)])
end

end
