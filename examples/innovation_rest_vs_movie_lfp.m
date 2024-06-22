addpath('../matlab/')
addpath('~/matlab/Violinplot-Matlab-master/'); % https://github.com/bastibe/Violinplot-Matlab/tree/master

% clear all

subjects = dir('/home/lparra/Dropbox/Intracranial/varx/version/varx_data/NS*');
condition = {'Despicable_Me_English_5min.mat', 'Resting_fixation.mat',};

% We will be plotting one subject at a time, and at the end add a summary.
% picking nice example subject to go last, so we have it in summayr figure.
last = 5;
subject_order = 1:length(subjects); subject_order(last)=[]; subject_order=[subject_order last];

for subj = 3 %subject_order

    % first load two conditions, and find electrodes where power 
    % does not change by more than 20dB, i.e. rule out contact problem or
    % movement artefacts
    for i=1:length(condition)
        load([subjects(subj).folder '/' subjects(subj).name '/' condition{i}],'hfa');
        y2{i} = var(hfa);
    end
    goodchannels = abs(db(y2{1}./y2{2})/2)<20; 

    % now analyze the data in the good channels
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
    
    drawnow

end


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