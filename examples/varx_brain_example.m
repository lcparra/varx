%% Here is code to estimate the impulse response of motion/saccade/fixation
% to BHA allowing for a "recurrent" input using a VARX model. 

addpath('../matlab/')
% stim_features = {'optical_flow','fixations'};
% str = {'flow','fixation'};
stim_features = {'fixations'};
str = {'fixation'};

videos = {...
    'Despicable_Me_English.mat', ...
    'Despicable_Me_Hungarian.mat', ...
    'Monkey1_Rep_1.mat', ...
    'Monkey2_Rep_1.mat', ...
    'Monkey5_Rep_1.mat', ...
    'The_Present_Rep_1.mat',...
    'The_Present_Rep_2.mat'};

data_path = '../data/Example_patient/'; channels=[31:80];
if ~exist(data_path)
    disp('Download test data here: https://osf.io/vc25t/?view_only=0aa33fa3d1df4114a2d040b01397c9d1')
    disp('unzip and place it here ../data/Example_patient/')
    return
end

tic
for vid=1:length(videos)
    
    video = videos{vid};
    
    load([data_path 'BHA/' video],'envelope','fs');
    load([data_path 'Stimuli/' video],stim_features{:});
       
    x = [];
    for f=1:length(stim_features)
        x = [x eval(stim_features{f})];
    end

    y = envelope(:,channels); 

    % convert power into something more normal distributed.
    y = log(y./std(y)+1);
       
    % subtract mean to avoid dealing with it in the regression
    x = x-mean(x); y=y-mean(y);

    % normalize scale for output so that penalty affects them equally
    y = y./std(y);
    x = x./std(x);

    X{vid}=x; Y{vid}=y;

    lengths(vid)=length(y);
        
end
t=toc; disp(['time to load data: ' num2str(t)]); 

%%
L=fs/2; % length of the FIR model filter
na=4; nb=20; % length of the VARX model filters
base = basis(fs/2,nb,'normal');
varx_model = varx(Y,na,X,base);

% derive the corresponding impulse responses
H = varx_trf(varx_model.B,varx_model.A,L);

% now compute MIMO FIR model (mTRF)
Rxx=0; Rxy=0;
for i=1:length(X)
    [Rxx_,Rxy_] = myxcorr(X{i},Y{i},L);
    Rxx = Rxx+Rxx_; Rxy = Rxy+Rxy_;
end
mTRF = inv(Rxx)*Rxy;
xdim=size(X{1},2);
ydim=size(Y{1},2);
mTRF = permute(reshape(mTRF,L,xdim,ydim),[1 3 2]);

% now lets do just an AR model
var_model = varx(Y,na);


%% some displays
figure(1)
varx_display(varx_model,fs=fs);   
% exportgraphics(gcf,'../figures/VARX-brain-example.png','Resolution',300)


%% compute conventional "functional connectivity", i.e. Pearson correlation matrix
tic; Nrand=100;
cc=corrcoef(cat(1,Y{:}));
for n=2:Nrand+1
    shifts = floor(rand(1,size(Y{1},2))*sum(lengths));
    cc(:,:,n)=corrcoef(mycircshift(cat(1,Y{:}),shifts));
end
cc_pval = max(mean(repmat(abs(cc(:,:,1)),1,1,Nrand)<cc(:,:,2:end),3),1/Nrand);
t=toc; disp(['time to compute shuffle stats for corr coeff: ' num2str(t)]); 

%% compare VARX B with mTRF estimate
figure(2); clf; h=[];
for i=1:size(varx_model.B,3) 

    h(end+1)=subplot(3,xdim,i); 
    imagesc((1:size(mTRF,1))/fs,1:size(y,2),mTRF(:,:,i)'); 
    ylabel('brain channels'); xlabel('lag (seconds)')
    title('mTRF Impulse response')
    clim([-1 1]*max(abs(mTRF(:)))); colorbar

    h(end+1)=subplot(3,xdim,xdim+i);
    imagesc((1:size(H,1))/fs,1:size(y,2),H(:,:,i)');
    ylabel('brain channels'); xlabel('lag (seconds)')
    title('VARX Impulse response')
    clim([-1 1]*max(abs(H(:)))); colorbar

end
% sublabel(h,-5,-20);    
% exportgraphics(gcf,'../figures/VARX-brain-example_model-compare.png','Resolution',300)

%% compare VARX with VAR results
subplot(3,2,1); imagesc(log10( var_model.A_pval)); ylabel('brain channels'); xlabel('brain channels'); title('p-value VAR');  axis equal;axis tight;  caxis([-3 0]); xlabel(colorbar,'log10(p)')
subplot(3,2,2); imagesc(log10(varx_model.A_pval)); ylabel('brain channels'); xlabel('brain channels'); title('p-value VARX'); axis equal;axis tight;  caxis([-3 0]); xlabel(colorbar,'log10(p)')
subplot(3,2,3); imagesc((varx_model.A_pval<0.01)-(var_model.A_pval<0.01));  ylabel('brain channels'); xlabel('brain channels'); title('change VARX to VAR'); axis equal;axis tight; 
subplot(3,2,4); imagesc(log10(cc_pval)); ylabel('brain channels'); xlabel('brain channels'); title('p-value Corr. Coeff.'); axis equal;axis tight; caxis([-3 0]); xlabel(colorbar,'log10(p)') 
subplot(3,1,3); plot([varx_model.s2 - var_model.s2]); ylabel('power difference (VARX - VAR)'); xlabel('brain channels'); title('Change in noise power ($\sum e^2/T$)','interpreter','latex');

%% help functions
function [X]=mycircshift(X,shifts)
for i=1:size(X,2)
    X(:,i)=circshift(X(:,i),shifts(i));
end
end


