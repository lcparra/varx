% Testing of GC with data generation following X->Y and a third variable Z
% which can be a confounder Z (X<-Z->Y), or collider (X->Z<-Y), or entirely
% indpependent (X Z Y). We estimate the Power of the log-likelihood
% Chi-square test by testing X->Y and the False discovery rate by testing
% Y->X. We test the effect of conditioning on Z (as proposed by Geweke in
% 1984), either only includeing the past (delayed Z) or also the present
% (instant Z), or not conditioning at all (none). We test two different
% data generation model coded by Behtash Babadi or myself. Behtash
% implements an equation error model, where I implement an output error
% model.
%
% Here my concludions:
%
%    In case that there is a collider, conditioning on it with instant
%    delay can introduce strong collider bias -- consistent with what
%    Judeas Pearl teaches us. Conditioning on the history of the collider
%    is safe.
%
%    Flip side, if there is an instant common cause as confound,
%    conditioning on it improves power. So, in general, use instant
%    de-confounding only if you are sure about cause and effect of the
%    de-confounding variable.
%
%    For smaller N false discovery goes up quite a bit for both models,
%    i.e. the test finds a GC link where there is none. This is worse for
%    the output error data generation model. So, trust GC only if you are highly
%    powered and are confident that the data generating model matches the
%    VAR model (equation error).
%
% Conclusion: "Grange Causality" is conventional linear regression
% following Pearl's instructions. The variables involved just so happen to
% be history of X and Y.
%
% Granger causal  x->y is the same as Pearl's: X(t-1) -> Y(t) <- Y(t-1)
%
% The only difference is in the statistical test. If I follow Pearl's
% instructions for this graph I get a p-value for each delay using
% conventional linear modeling -- in matlab, I use fitlm() for this.
% Whereas with the likelihood ratio test of GC, here we do a single test
% for the entire history, probably gaining in power. The role and treatment
% of confounders and colliders is exactly the same. So really, to me now GC
% is just a causal linear model on the history. The test is correct if the
% graph is correct. Period. You get direction only because of the
% assumption that the present depends on the past. When you have no strong
% confidence about how instant effects operate, I will stay away from GC to
% make statements about direction of an effect.

Q=3; % delays in model and estimation
N=5000; % number of samples

Nrand=1000; % number of trials to estiamte power and false discovery rate

% things to try
confound = {'collider', 'independent','common cause'};
control = {'none','instant','delayed'};
model = {'Equation Error Model','Output Error Model'};

for time_reversed=0:1

    figure(time_reversed+1); clf

    Power = zeros(length(confound),length(control),length(model));
    FDrate= zeros(length(confound),length(control),length(model));
    for imodel=1:length(model)
        for icont=1:length(control)
            for iconf=1:length(confound)

                true_discovery = zeros(Nrand,1);
                false_discovery = zeros(Nrand,1);
                for irand=1:Nrand % test how many times I get a false alarm
                    switch model{imodel}
                        case 'Equation Error Model' % Behtash Babadi
                            % coefficients here starting with daly 1
                            byx = 1;    % effect of X on Y
                            ax = 0.1*rand(Q,1); % AR coefficient of outcome X
                            ay = 0.1*rand(Q,1); % AR coefficient of outcome Y
                            x = randn(N,1); % generate some random "cause" X
                            y = zeros(N,1);
                            for k=Q+1:N
                                y_e = 0; % effective history of Y at time k
                                x_e = 0; % effective hisotyr of X at time k
                                for l=1:Q
                                    y_e = y_e + ay(l)*y(k-l);
                                    x_e = x_e + ax(l)*x(k-l);
                                end
                                y(k)= y_e + byx*x_e  + randn(1,1); % adding y_e and x_e and some obseevation noise
                            end
                        case 'Output Error Model' % Lucas Parra'
                            % coefficients here starting with daly 0, as in filter() function
                            byx =[0; rand(Q,1)]; % effect of X on Y, can't have instant effect, else can act as confound
                            ax = [1; 0.1*rand(Q,1)]; % AR coefficient of outcome X
                            ay = [1; 0.1*rand(Q,1)]; % AR coefficient of outcome Y
                            x = filter(  1,ax,  randn(N,1)); % generate some random "cause" X, plus confound
                            y = filter(byx,ay,x)+randn(N,1); % add it to some AR process, plus confound
                    end
                    switch confound{iconf}
                        case 'common cause'
                            z = randn(N,1); % comon cause confounder
                            y = y + z; x = x + z; % add instant effect confound
                        case 'collider'
                            z = x + y + 0.1*randn(N,1); % simulate collider
                        case 'independent'
                            z = randn(N,1); % neither cause or effect to x, y
                    end

                    if time_reversed
                        x=flipud(x);y=flipud(y);z=flipud(z);
                    end

                    true_discovery(irand)  = granger_test(x,y,Q,z,control{icont});
                    false_discovery(irand) = granger_test(y,x,Q,z,control{icont});

                end

                Power (iconf,icont,imodel) = mean(true_discovery <0.05);
                FDrate(iconf,icont,imodel) = mean(false_discovery<0.05);

            end
        end

        subplot(2,2,(imodel-1)*2+1);
        bar(Power (:,:,imodel)); set(gca,'xticklabel',confound);
        title(model{imodel}); xlabel('Confound'); ylabel('Power');
        subplot(2,2,(imodel-1)*2+2);
        bar(FDrate(:,:,imodel)); set(gca,'xticklabel',confound); hold on
        title(model{imodel}); xlabel('Confound'); ylabel('False Discovery Rate');
        ax = axis; plot(ax(1:2),[0.05 0.05],'--'); hold off
        legend([control '\alpha=0.05'])
        drawnow
    end

    if time_reversed
        suptitle('Granger estimate with time reversal')
        saveas(gcf,'figures/granger_simulation_timereversed.png')
    else
        saveas(gcf,'figures/granger_simulation.png')
    end
end

