
%% VARX identification on real data
clear all
load SteamEng % arbitrary data to demo the code, not to test performance!
x=[Pressure,MagVolt]; % some multidimentional input
y=[GenVolt,Speed];
na=10; nb=20; gamma=0.3;

[Aest,Best,A_pval,B_pval] = varx(y,na,x,nb,gamma);
A_pval,B_pval

yest = varx_simulate(Best,Aest,x,y);

figure(3)
show_results(x,y,yest,Aest,Best,A_pval,B_pval);
sgtitle('Real data exampe')


%% generate a data with known VARX model and compare varx result with matlab's varm estimate
clear all

% define VARX model. Where there are all zeros, that path does not
% contribute
A(:,:,1) = [[0.3 -0.5 0]',[0 0 0]']; A(:,:,2) = [[-.5 0.4 0]',[0.5 -.7 0]'];
B = [[2 0]',[0 0]']; % single input. the varm estimate only works for zero delay
[nb,ydim,xdim] = size(B);
[na,ydim,ydim] = size(A);

% simulate 
figure(1); clf
gamma = 0;
T = 1000;
x = randn(T,xdim)+2;
[y,e] = varx_simulate(B,A,x,1); y=y+2; 

% estimate VARX model
[Aest,Best,A_pval,B_pval] = varx(y,na,x,nb,gamma);

[yest,eest] = varx_simulate(Best,Aest,x,y);

R2=1-sum(eest.^2)./sum(y.^2)

% estimate model using Econometric toolbox' mvar
ydim = size(A,2);
Mdl = varm(ydim,na);
EstMdl = estimate(Mdl,y,'X',x);

% some displays
figure(1)
show_results(x,y,yest,Aest,Best,A_pval,B_pval);
sgtitle('Simulated VARX and estimation')

% compare estimate to truth
figure(2); clf
for t=1:na, for i=1:ydim, for j=1:ydim, AR{i,j}(t) = EstMdl.AR{t}(i,j); end; end; end;
plot(A(:),Aest(:),'.'); hold on;
plot(A(:),stack([AR{:}]),'.');
plot(B(:),Best(:),'x')
plot(B(1,:),EstMdl.Beta(:),'x')
plot([-1 1],[-1 1]); hold off
ylim([-1.5 1.5]); %axis equal; axis tight
legend('our AR estimate','matlab AR estimte','our MA estimate','matlab MA estimte','Location','NorthWest')
xlabel('true value'); ylabel('estimate')

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
            gamma = lambda(i)/sqrt(T(nt)-na*ydim-nb*xdim);
            [Aest,Best,A_pval,B_pval] = varx(y(1:T(nt),:),na,x(1:T(nt),:),nb,gamma);
            Apresent(i,nt,nrand)=A_pval(1,2); 
            Amissing(i,nt,nrand)=A_pval(2,1);
            Bpresent(i,nt,nrand)=B_pval(2,2); 
            Bmissing(i,nt,nrand)=B_pval(1,2);
            [yest,eest] = varx_simulate(Best,Aest,x,y);
            E_train(i,nt,nrand)= mean(sum(eest(1:T(nt)    ,:).^2)./sum(y(1:T(nt)    ,:).^2));
            E_test (i,nt,nrand)= mean(sum(eest(T(nt)+1:2*T(nt),:).^2)./sum(y(T(nt)+1:2*T(nt),:).^2));
        end
    end
    subplot(2,2,1)
    errorbar(lambda, mean(E_train,3),std(E_train,[],3)/sqrt(size(E_train,3))); 
    ylabel('Training Error'); xlabel('\lambda')
    subplot(2,2,2)
    errorbar(lambda, mean(E_test, 3),std(E_test ,[],3)/sqrt(size(E_test ,3))); 
    ylabel('Test Error'); 
    clear str; for i=1:length(T), str{i}=['T=' num2str(T(i))]; end; 
    legend({str{:},'location','northeast'})
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

%saveas(gcf,'figures/effect_of_regularization.png')


%% try repeating this 1000 times and see if the pvalues are correct
clear all
A(:,:,1) = [[0.9 -0.5 0]',[0 0 0]']; A(:,:,2) = [[-.5 0.4 0]',[0.5 -.7 0]'];
B = [[1 0 0 0 -1 0 0 0 0 1 0 0 -1 ]',[0 0 0 0 0 0 0 0 0 0 0 0 0 ]']; % the varm estimate only works for zero delay
[nb,ydim,xdim] = size(B);
[na,ydim,ydim] = size(A);
clear A_pval B_pval
for n=1000:-1:1
    x = randn(100,size(B,3)); 
    y = varx_simulate(B,A,x,1); 
    x(30,1) = NaN; % test if adding NaN works
    [~,~,A_pval(:,:,n),B_pval(:,:,n)] = varx(y,na,x,nb,0.3);
end
mean(A_pval<0.05,3)
mean(B_pval<0.05,3)

%% now do it with some larger models
clear all
na = 2; ydim=3;
nb = 2; xdim=1;
for i=1:ydim, 
    for j=1:ydim, A(:,i,j) = 0.1*randn(na,1); end; 
    for j=1:xdim, B(:,i,j) = randn(nb,1); end; 
end
% set two of the paths to zero, to see if we correctly identify them at alpha<0.05
A(:,2,end)=0;
B(:,1,end)=0;
clear A_pval B_pval
for n=100:-1:1
    x = randn(2000,xdim);
    y = varx_simulate(B,A,x,1);
    [~,~,A_pval(:,:,n),B_pval(:,:,n)] = varx(y,na,x,nb,0.1);
end
mean(A_pval<0.05,3)
mean(B_pval<0.05,3)'


%% now try using basis functions
clear A B
A(:,:,1) = [[0.9 -0.5 0]',[0 0 0]']; A(:,:,2) = [[-.5 0.4 0]',[0.5 -.7 0]'];
B = [[1 0.5 0 -0.5]',[0 0 0 0]']; % basis coefficients
[nb,ydim,xdim] = size(B);
[na,ydim,ydim] = size(A);
base = basis(30,nb,'normal');
gamma = 1;

x = randn(100,xdim);
y = varx_simulate(tensorprod(base,B,2,1),A,x,1);
[Aest,Best,A_pval,B_pval] = varx(y,na,x,base,gamma)

yest = varx_simulate(tensorprod(base,Best,2,1),A,x,y);
% x(30,1) = NaN; % test if adding NaN works

figure(4)
show_results(x,y,yest,Aest,Best,A_pval,B_pval);

figure(5); clf
plot(A(:),Aest(:),'o'); hold on
plot(B(:),Best(:),'o'); 
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

for n=5000:-1:1
    x = randn(300,xdim);
    % generate using basis function
    y = varx_simulate(tensorprod(base,B,2,1),A,x,1);
    % estimate using basis functions
    [~,Best_b,A_pval_b(:,:,n),B_pval_b(:,:,n)]  = varx(y,na,x,base);    
    Best_b = tensorprod(base,Best_b,2,1);
    % estimate all parameters
    [~,Best,A_pval(:,:,n),B_pval(:,:,n)]  = varx(y,na,x,size(base,1));
end
%%
clf
haxis(1)=subplot(1,3,1);
FDrate(1,1) = mean(A_pval(2,1,:)<0.05,3);
FDrate(1,2) = mean(B_pval(2,1,:)<0.05,3);
FDrate(2,1) = mean(A_pval_b(2,1,:)<0.05,3);
FDrate(2,2) = mean(B_pval_b(2,1,:)<0.05,3);
bar(FDrate); hold on
ylabel('False discovery rate'); ylim([0 0.07])
ax=axis; plot(ax(1:2),[0.05 0.05],':k'); hold off
legend('AR link','MA link','\alpha=0.05','location','northwest')
set(gca,'XTickLabel',{'free','basis'})
xlabel('Estimate')
haxis(2)=subplot(1,3,2);
for i=1:2
    plot(Best(:,i),'b'); hold on;
    plot(Best_b(:,i),'r');
    plot(tensorprod(base,B(:,i),2,1),'g'); 
end
legend('free estimate','basis estimate','true MA filter')
hold off
xlabel('lag (samples)')
ylabel('MA filter B')
haxis(3)=subplot(1,3,3); plot(base); 
xlabel('lag (samples)'); ylabel('basis functions')
sublabel(haxis,-5,-40);

saveas(gcf,'figures/effect_of_basis.png')

%% demo modeling just the AR part
clear all
A(:,:,1) = [[0 0 0]',[0. -0.5 0]']; A(:,:,2) = [[-.5 0.4 0]',[0.5 -.7 0]'];
ydim = size(A,2); xdim=0;
B=zeros(0,ydim,0);
na=size(A,1); nb=size(B,1); gamma=0.5;
x=randn(200,xdim); % no external input
[y,e] = varx_simulate(B,A,x,1);

[Aest,Best,A_pval,B_pval] = varx(y,na);
A_pval

[yest,e] = varx_simulate(Best,Aest,x,y);

% estimate the Ainv
[~,Ainv] = varx_trf(Best,Aest,100);

figure(4)
% abusing show_result function to also show Ainv instead of Best
show_results(x,y,yest,Aest,Ainv,A_pval,A_pval);
% correcting the abuse by correcting titles
for i=1:ydim, subplot(4,ydim,3*ydim+i); title(['AR inverse ' num2str(i)]); end

figure(5); clf
plot(A(:),Aest(:),'o'); hold on
plot(B(:),Best(:),'o'); 
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


A_pval = nan(ydim,ydim); % off dianoal elements are not computed, make sure that is clear. 
for i=1:ydim
    [Aest(:,i,i),Best(:,i,:),A_pval(i,i),B_pval(i,:)] = varx(y(:,i),na,x,nb);
end
A_pval,B_pval

[yest,e] = varx_simulate(Best,Aest,x,y);

figure(4)
show_results(x,y,yest,Aest,Best,A_pval,B_pval);

figure(5); clf
plot(A(:),Aest(:),'o'); hold on
plot(B(:),Best(:),'o'); 
plot([-1 1],[-1 1],':k'); hold off
xlabel('True values'); ylabel('Estimated values')
legend('MA coeeficients','location','northwest')

% derive the corresponding impulse responses
figure(6); clf
H = varx_trf(Best,Aest,20);
imagesc(H(:,:))

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
[Aest,~,A_pval_missing_input,B_pval] = varx(y,na); A_pval_missing_input

% simulate varx without the external imput (I know, its cumbersome) 
yest = varx_simulate(zeros(0,ydim,0),Aest,zeros(size(y,1),0),y);

% some displays
figure(1)
show_results([],y,yest,Aest,[],A_pval_missing_input,[]);
sgtitle('Simulated VARX and estimation')

% compare estimate to truth
figure(2); clf
plot(A(:),Aest(:),'.'); hold on;
plot([-1 1],[-1 1]); hold off
xlabel('true value'); ylabel('estimate')
sgtitle('Simulated VARX and estimation')

% estimate VARX model now with the input included
[Aest,Best,A_pval_input_included,B_pval] = varx(y,na,x,nb); A_pval_input_included

yest = varx_simulate(Best,Aest,x,y);

% some displays
figure(3)
show_results(x,y,yest,Aest,Best,A_pval_input_included,B_pval);
sgtitle('Simulated VARX and estimation')

% compare estimate to truth
figure(4); clf
plot(A(:),Aest(:),'.'); hold on;
plot([-1 1],[-1 1]); hold off
xlabel('true value'); ylabel('estimate')
sgtitle('Simulated VARX and estimation')


%% model long smooth AR coefficients, see if we recover them correctly. 
% make smooth A
clear A
A(:,:,1) = [[0.1 -0.5 0]',[.1 0 0]']; A(:,:,2) = [[-0.1 0 0]',[0.05 -.7 0]'];
[na,ydim,ydim] = size(A);
base = basis(30,na,'normal');
A = tensorprod(base,A,2,1);
[na,ydim,ydim] = size(A)

% make dummy input B, x
xdim=0;
B=zeros(0,ydim,0);
na=size(A,1); nb=size(B,1);
x=zeros(2000,xdim); % no external input

% simulate 
[y,e] = varx_simulate(B,A,x,1);

% estimate VARX model 
gamma = 0.0
[Aest,~,A_pval,~] = varx(y,na,[],[],gamma); A_pval

% simulate equaiton error model with the estimated Aest
yest = varx_simulate(B,Aest,x,y);

% some displays
figure(1)
show_results([],y,yest,Aest,[],A_pval,[]);
sgtitle('Simulated VARX and estimation')
plot(A(:,:)); hold on
plot(Aest(:,:)); hold off

% compare estimate to truth
figure(2); clf
plot(A(:),Aest(:),'.'); hold on;
plot([-1 1],[-1 1]); hold off
xlabel('true value'); ylabel('estimate')
sgtitle('Simulated VARX and estimation')




% ----------------- result display function --------------------------
function show_results(x,y,yest,A,B,A_pval,B_pval)

clf
ydim = size(A,2);
xdim = size(B,3);
for i=1:ydim
    subplot(4,2,1); plot(x); title('External Inputs');    
    subplot(4,2,2); plot(y); title('Recursive Input');    
    subplot(4,ydim,ydim*1+i); 
    plot([y(:,i) yest(:,i)]) % compare original to estimated output
    title(['Target and estimated output ' num2str(i)])
    subplot(4,ydim,ydim*2+i); 
    for j=1:ydim % show AR filters
        plot(A(:,i,j),'LineWidth',0.5+(A_pval(i,j)<0.05)); hold on; 
    end
    hold off; title(['IIR filter estimates ' num2str(i)])
    subplot(4,ydim,ydim*3+i); 
    for j=1:xdim % show MA filters
        plot(B(:,i,j),'LineWidth',0.5+(B_pval(i,j)<0.05)); hold on; 
    end
    hold off; title(['FIR filter estimates ' num2str(i)])
    xlabel('time delay (samples)')
end
end

