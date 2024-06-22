addpath('../matlab/')
addpath('~/matlab/Violinplot-Matlab-master/'); % https://github.com/bastibe/Violinplot-Matlab/tree/master

% clear all

subjects = dir('/home/lparra/Dropbox/Intracranial/varx/version/varx_data/NS*');
condition = {'Despicable_Me_English_5min.mat', 'Resting_fixation.mat',};

% We will be plotting one subject at a time, and at the end add a summary.
% picking nice example subject to go last, so we have it in summayr figure.
last = 4;
subject_order = 1:length(subjects); subject_order(last)=[]; subject_order=[subject_order last];

for subj = subject_order

    % first load two conditions, and find electrodes where power 
    % does not change by more than 20dB, i.e. rule out contact problem or
    % movement artefacts
    for i=1:length(condition)
        load([subjects(subj).folder '/' subjects(subj).name '/' condition{i}],'hfa');
        y2{i} = var(hfa);
    end
    goodchannels = abs(db(y2{1}./y2{2})/2)<20; 

    % now analyze the data in the good channels
    figure(1); clf; tiledlayout(3,3)
    for i=1:length(condition)
    
        warning off % i know some of these variables are missing in rest condition
        load([subjects(subj).folder '/' subjects(subj).name '/' condition{i}],'lfp','fs','fixations','film_cuts','audio_env');
        warning on

        % The order of loading is important, first movie then rest, so that
        % film_cuts and audio_env are defined also for rest. We want to
        % regress out the same signal in both conditions to have same
        % number of free parameters in both conditions.
        x = [fixations film_cuts audio_env]; 
        y = lfp(:,goodchannels); clear lfp;
        
        % subtract mean to avoid dealing with it in the regression
        x = x-mean(x); y=y-mean(y);

        y2{i} = var(y);
        [Py{i},fbin]=pwelch(y,fs,[],[],fs);

        % normalize scale for output so that penalty affects them equally
        y = y./sqrt(y2{1});
        x = x./std(x);

        na = 4; nb=fs/2; lambda=0;
        model{i} = varx(y,na,x,nb,lambda);

        [~,e] = varx_simulate(model{i}.B,model{i}.A,x,y);
        [Pe{i},fbin]=pwelch(e,fs,[],[],fs);

        h1(i)=nexttile((i-1)*3+1,[1 2]); 
        plot(model{i}.B_Rvalue); axis tight
        ylabel('Effect size R')
        if i==1, title('movie'); else, title('rest'); end;

    end

    % remember for each subject change in power for signal and noise
    y2_diff{subj} = db(y2{1})-db(y2{2});
    s2_diff{subj} = db(model{1}.s2)-db(model{2}.s2);

    % display change in signal and noise power resolved by frequency for all electrodes
    figure(2), tiledlayout(2,3), 
    nexttile(1,[1 2]); 
    imagesc(1:size(Py{1},2),fbin,db(abs(Py{1}))-db(abs(Py{2}))); h(1)=gca; pos=h(1).Position;
    colorbar('westoutside'); clim([-10 10]); axis xy; h(1).Position=pos;
    ylabel('Frequency (Hz)'); title('y^2 Power spectrum change');  
    nexttile(4,[1 2]), 
    imagesc(1:size(Pe{1},2),fbin,db(abs(Pe{1}))-db(abs(Pe{2}))); h(2)=gca; pos=h(2).Position;
    colorbar('westoutside'); clim([-10 10]); axis xy; h(2).Position=pos;
    ylabel('Frequency (Hz)'); title('e^2 Power spectrum change'); xlabel('Channels')
    

    % report change in power for responsive and non-responsive electrodes
    responsive=sum(model{1}.B_pval<0.01/size(y,2)/size(x,2),2)>0; % bonferony corrected
    Change_responsive    {subj} = noise_diff( responsive);
    Change_not_responsive{subj} = noise_diff(~responsive);
    if sum(responsive)>0
        pval = ranksum(noise_diff(responsive),noise_diff(~responsive));
    else pval=1; end
    display([subjects(subj).name ...
        ': responsive '   num2str(median(noise_diff( responsive)),2), ...
        ' non responsive ' num2str(median(noise_diff(~responsive)),2), ...
        ', p=' num2str(pval,2) ...
        ' N=' num2str(sum(responsive)) ' responsive electrodes.'])
    drawnow


end

%% summary of change in power of relative noise in responsive and non responsive electrodes
figure(1); 
h1(4)=nexttile(3,[3 1])
N=length(Change_responsive);
clear tmp1; for i=1:N, tmp1(i)=median(Change_responsive{i}); end; tmp1=tmp1';
clear tmp2; for i=1:N, tmp2(i)=median(Change_not_responsive{i}); end; tmp2=tmp2'; 
pval = ranksum(tmp1,tmp2);
display(['Change in power ratio responsive vs non responsive channesl, p=' num2str(pval,2) ', N=' num2str(N) ' subjects.'])

violinplot([tmp1;tmp2],[repmat("responsive",size(tmp1,1),1);repmat("not responsive",size(tmp2,1),1)]); 
ax=axis; plot(ax(1:2),[0 0],'k'); hold off 
title('Median Channel')
ylabel('e^2/y^2 change (dB)')
sublabel(h1,-15,-30);
saveas(gca,'../figures/noise-power-relative-change_responsive.png')

%% show change in power of signal and noise in all electrodes
figure(2); 
N=length(y2_diff);
nexttile(3,[2 1]);  
clear tmp1; for i=1:N, tmp1(i)=median(y2_diff{i}); end; tmp1=tmp1';
clear tmp2; for i=1:N, tmp2(i)=median(s2_diff{i}); end; tmp2=tmp2'; 
violinplot([tmp1;tmp2],[repmat("y^2",size(tmp1,1),1);repmat("e^2",size(tmp2,1),1)]); 
ax=axis; plot(ax(1:2),[0 0],'k'); hold off 
ylim([-1 1]*max(abs(ax(3:4))))
title('Median Channel'); ylabel('Power change (dB)')
display(['Change in y^2, p=' num2str(signrank(tmp1)) ' Subject N=' num2str(N)])
display(['Change in e^2, p=' num2str(signrank(tmp2)) ' Subject N=' num2str(N)])
h(3)=gca;
sublabel(h,-15,-30);

saveas(gca,'../figures/noise-power-change.png')