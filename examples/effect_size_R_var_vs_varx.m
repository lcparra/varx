% Here we want to test what happes to effect size R of recurrent
% connections when we estimate the data ignoring an input (VAR),
% including the input (VARX) or using the same number of free parameters
% but using a random input that does not match the data generation (VARXn).

% Conclusion: VAR and VARXn behave the same. Adding the correct X (i.e.
% VAR->VARX) modulates R for in-connections to the X-responsive node,
% either increasing of decreasing it according to the sign of that
% connection.

addpath('../matlab/')
addpath('~/matlab/Violinplot-Matlab-master/'); % https://github.com/bastibe/Violinplot-Matlab/tree/master


%% first check on a simulated VARX toy model
A(:,1,:) = [[0.3 -0.5 0]',[0 0 0    ]', [-0.5 0 0]']; 
A(:,2,:) = [[.5 0.4  0]',[0.5 -.7 0]', [0 0 0]'];
A(:,3,:) = [[-1  0   0]',[0 0 0]', [0.5 0 0]'];
B = [[1 0]',[0 0]',[0 0]']; % single input.
[nb,ydim,xdim] = size(B);
[na,ydim,ydim] = size(A);

% simulate VARX with single input
T = 1000;
x = randn(T,xdim);
[y,e] = varx_simulate(B,A,x,1);  

% estimate VARX model
model(1,1) = varx(y,na);
model(2,1) = varx(y,na,x,nb);
model(3,1) = varx(y,na,randn(T,xdim),nb);

% simulare VAR with no input
[y,e] = varx_simulate(zeros(0,ydim,0),A,zeros(T,0),1);
model(1,2) = varx(y,na);
model(2,2) = varx(y,na,x,nb); % regress out even if not the correct input, just so the DF are the same

figure(1)
show_results(1,model,1,'Toy example')



%% repeat this for LFP data

subjects = dir('/home/lparra/Dropbox/Intracranial/varx/version/varx_data/NS*');
condition = {'Despicable_Me_English_5min.mat', 'Resting_fixation.mat',};

% We will be plotting one subject at a time, and at the end add a summary.
% picking nice example subject to go last, so we have it in summayr figure.
last = 5;
subject_order = 1:length(subjects); subject_order(last)=[]; subject_order=[subject_order last];

for subj =  1 % subject_order

    % now analyze the data in the good channels
    for i=1:length(condition)
    
        warning off % i know some of these variables are missing in rest condition
        load([subjects(subj).folder '/' subjects(subj).name '/' condition{i}],'lfp','fs','fixations','film_cuts','audio_env');
        warning off

        % The order of loading is important, first movie then rest, so that
        % film_cuts and audio_env are defined also for rest. We want to
        % regress out the same signal in both conditions to have same
        % number of free parameters in both conditions.
        x = [fixations film_cuts audio_env]; 
        y = lfp(:,:); clear lfp;
        
        % subtract mean to avoid dealing with it in the regression
        x = x-mean(x); y=y-mean(y);
        y2{i} = var(y);

        % normalize scale for output so that penalty affects them equally
        y = y./sqrt(y2{1});
        x = x./std(x);

        na = 4; nb=fs/2;
       
        % estimate  models
        model(1,i,subj) = varx(y,na);
        model(2,i,subj) = varx(y,na,x,nb);
        model(3,i,subj) = varx(y,na,[x(:,1) randn(size(x,1),2)],nb);

    end
    
    show_results(2,model(:,:,subj),0.1,'Movie vs Rest')
    
end
return

%%
figure(2)
clf
for subj=subject_order
    mdiff(subj,1,:)=median_diff_sig_R(model(2,1,subj),model(2,2,subj)); 
    mdiff(subj,2,:)=median_diff_sig_R(model(2,1,subj),model(3,1,subj));
    mdiff(subj,3,:)=median_diff_sig_R(model(3,1,subj),model(2,2,subj));
end
for i=1:2
subplot(2,3,3*(i-1)+1); violinplot(mdiff(:,1,i));[h,p]=ttest(mdiff(:,1,i)); xlabel('VARX Stim - VARX noStim');  title(['p=' num2str(p,2)]); ax=axis;plot(ax(1:2),[0 0],'k')
subplot(2,3,3*(i-1)+2); violinplot(mdiff(:,2,i));[h,p]=ttest(mdiff(:,2,i)); xlabel('VARX Stim - VARXn Stim');    title(['p=' num2str(p,2)]); ax=axis;plot(ax(1:2),[0 0],'k')
subplot(2,3,3*(i-1)+3); violinplot(mdiff(:,3,i));[h,p]=ttest(mdiff(:,3,i)); xlabel('VARXn Stim - VARXn noStim'); title(['p=' num2str(p,2)]); ax=axis;plot(ax(1:2),[0 0],'k')
end
subplot(2,3,1); ylabel('R value Median diff of signifficant channels')
subplot(2,3,4); ylabel('Fraction Median diff of signifficant channels')
saveas(gcf,'../figures/VARX Movie vs VAR Rest vs VAR Movie.png')

function show_results(fig,model,cscale,name)
  figure(fig); clf
  subplot(4,3,1); imagesc(model(1,1).A_Rvalue); clim(cscale*[0 1]); title('VAR: R A')
  subplot(4,3,4); imagesc(model(2,1).A_Rvalue); clim(cscale*[0 1]); title('VARX: R A')
  subplot(4,3,7); imagesc(model(3,1).A_Rvalue); clim(cscale*[0 1]); title('VARXn: R A')
  subplot(4,3,2); imagesc(model(1,1).B_Rvalue); clim(cscale*[0 1]); title('VAR: R B')
  subplot(4,3,5); imagesc(model(2,1).B_Rvalue); clim(cscale*[0 1]); title('VARX: R B')
  subplot(4,3,8); imagesc(model(3,1).B_Rvalue); clim(cscale*[0 1]); title('VARXn: R B')
  subplot(4,3,3); imagesc(model(1,1).A_Rvalue-model(3).A_Rvalue); clim(cscale*0.1*[-1 1]); title('VAR -VARXn')
  subplot(4,3,6); imagesc(model(2,1).A_Rvalue-model(1).A_Rvalue); clim(cscale*0.1*[-1 1]); title('VARX -VAR')
  subplot(4,3,9); imagesc(model(2,1).A_Rvalue-model(3).A_Rvalue); clim(cscale*0.1*[-1 1]); title('VARX-VARXn')
  mdiffR=median_diff_sig_R(model(2,1),model(3,1)); 
  disp([name ': Stim median R VARX - VARXn: ' num2str(mdiffR,2)])

  subplot(4,3,10); imagesc(model(2,1).A_Rvalue); clim(cscale*[0 1]); title('VARX Stim')
  subplot(4,3,11); imagesc(model(2,2).A_Rvalue); clim(cscale*[0 1]); title('VAR noStim')
  subplot(4,3,12); imagesc(model(2,1).A_Rvalue-model(1,2).A_Rvalue); clim(cscale*0.1*[-1 1]); title('VARX Stim - VAR noStim')
  mdiffR=median_diff_sig_R(model(2,1),model(2,2)); 
  disp([name ': median R VARX Stim - VARXn noStim: ' num2str(mdiffR,2)])
  mdiffR=median_diff_sig_R(model(3,1),model(2,2));
  disp([name ': median R VARXn Stim - VARXn noStim: ' num2str(mdiffR,2)])

  suptitle(name)
  
  drawnow
end

function y = median_diff_sig_R(model1,model2)

  Rdiff = model1.A_Rvalue-model2.A_Rvalue;
  mask = model1.A_pval<0.0001 | model2.A_pval<0.0001;
  y(1) = median(Rdiff(find(mask.*(ones(size(Rdiff))-eye(size(Rdiff))))));
  y(2) = (sum(model1.A_pval(:)<0.0001)-sum(model2.A_pval(:)<0.0001))/sum(mask(:));

end

