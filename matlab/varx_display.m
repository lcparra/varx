function [Graph,G_plot]=varx_display(m,o)
% [Graph,G_plot]=varx_display(model,options)
% Displays the model effect sizes and model filter parameters. Example:
% 
% model = varx(y,na,x,nb);
% varx_display(model,xname={},yname={},duration=0,fs=1,threshold=0.001);
% 
% The argument model is a structure with the following stucture elements
% A,B,A_pval,B_pval,A_Deviance,B_Deviance,T.base as created by model =
% varx(). 
% 
% options are name=value optional arguments:  yname,xname,duration,fs,threshold
%
% yname,xname are cell arrays that specify the names of input x and output
% y in the model. If they are provided, then model structure is displayed
% as a circular graph, and Graph and G_plot can be used to plot with:
%
% plot(Graph,'LineWidth',G_plot.width,'XData',G_plot.xdata,'YData',G_plot.ydata,...
%            'EdgeColor',G_plot.color,'NodeColor',G_plot.nodecolor);
%
% If variable ynames is not provided, everything is shown as matrix
% instead. Suggest to use this when there are a lot of variables in the
% model.
%
% duration is a scalar indicating for how long to compute (display) the
% impulse respones. IR are from every variable to every other variable. If
% it is omited then the IR is not displayed. 
%
% fs is the sampling rate. If given, then filters and impulse response are
% displayed in seconds.
%
% threshold is the cutoff pvalue to display links. Only links (effect
% size) below this cutoff will be displayed. The effect size is the R-value
% (see varx.m for detail).

% March 11, 2024, Lucas Parra
% March 16, 2024, changed to display R as effect size. added threshold as argument, and added options arguments 

% define arguments with defaults for the options variables o
arguments
    m 
    o.yname cell = {}     % show as matrix 
    o.xname cell = {}     
    o.duration double = 0 % no impulse response is shown
    o.fs double = 1       % lags displayed in samples. 
    o.threshold double = 0.001 % only links with lower values shown
end

clf

[nb,~,Dx] = size(m.B);
[na,~,Dy] = size(m.A);

if ~isempty(o.yname)

    pval = [m.A_pval m.B_pval; ones(Dx,Dy+Dx)];

    % sqrt of Coeficient of variation, i.e. R-square=1-exp(-D/T)
    Effect = [m.A_Deviance m.B_Deviance; zeros(Dx,Dy+Dx)]/m.T;
    Rvalue = sqrt(1-exp(-Effect)); 

    subplot(2,1,1)

    % plot connectivity as a graph
    node_name = [o.yname(:)' o.xname(1:Dx)];

    % set graph visuals
    h=plot(digraph(zeros(Dy)),'Layout','circle');
    G_plot.ydata=[get(h,'YData') linspace(.8,-0.8,Dx)];
    G_plot.xdata=[get(h,'XData') -2*ones(1,Dx)];
    G_plot.nodecolor = [repmat([0 0 0.75],Dy,1); repmat([0.75 0 0],Dx,1)];

    % generate the graph and set more graph visuals
    [s,t]=meshgrid(1:Dy+Dx,1:Dy+Dx);
    Adj = pval<o.threshold;
    Graph = digraph(s(Adj),t(Adj),Rvalue(Adj),node_name,'omitselfloops');
    G_plot.width = 50*Graph.Edges.Weight;
    if numedges(Graph)
        for i=1:numedges(Graph)
            if any(strcmp(Graph.Edges.EndNodes{i,1},o.xname))
                G_plot.color(i,:) = [1 0.25 0.25];
            else
                G_plot.color(i,:) = [0.25 0.25 1];
            end
        end
    else
        G_plot.color = 'none';
    end
    % now plot the Graph
    plot(Graph,'LineWidth',G_plot.width,'XData',G_plot.xdata,'YData',G_plot.ydata,'EdgeColor',G_plot.color,'Nodecolor',G_plot.nodecolor);
    axis equal; axis off;

    % make and plot the TRFs
    for i=1:Dy
        subplot(4,Dy+Dx,2*(Dy+Dx)+i);
        plot((1:na)/o.fs,m.A(:,:,i).*shiftdim(pval(1:Dy,i)<0.0001,1)); stack
        set(gca,'ytick',[]); title(node_name{i})
        if o.fs==1, xlabel('samples'); else, xlabel('seconds'); end
    end
    for i=1:Dx
        subplot(4,Dy+Dx,2*(Dy+Dx)+Dy+i);
        plot((1:nb)/o.fs*1,m.B(:,:,i).*shiftdim(pval(1:Dy,Dy+i)<0.0001,1));
        set(gca,'ytick',[]); title(node_name{Dy+i}); stack
        if o.fs==1, xlabel('samples'); else, xlabel('seconds'); end
    end
    subplot(4,Dy+Dx,2*(Dy+Dx)+1); axis on; ylabel('AR,MA filters')
    pos=get(gca,'Position');
    h=legend(node_name{1:Dy},'Location','westoutside');
    set(gca,'Position',pos);
    title(h,'Effect on')

    % show the impulse responses for duration
    if o.duration
        [H,Ainv] = varx_trf(m.B,m.A,o.duration*o.fs);
        for i=1:Dy
            subplot(4,Dy+Dx,3*(Dy+Dx)+i); plot((1:o.duration*o.fs)/o.fs,Ainv(:,:,i));
            set(gca,'ytick',[]); title(node_name{i}); stack
            if o.fs==1, xlabel('samples'); else, xlabel('seconds'); end
        end
        for i=1:Dx
            subplot(4,Dy+Dx,3*(Dy+Dx)+Dy+i); plot((1:o.duration*o.fs)/o.fs,H(:,:,i));
            set(gca,'ytick',[]); title(node_name{Dy+i}); stack
            if o.fs==1, xlabel('samples'); else, xlabel('seconds'); end
        end
        subplot(4,Dy+Dx,3*(Dy+Dx)+1); axis on; ylabel('Impulse response')
        pos=get(gca,'Position');
        h=legend(node_name{1:Dy},'Location','westoutside');
        set(gca,'Position',pos);
        title(h,'Response of')
    end

else

    % sqrt of Coeficient of variation, i.e. R-square
    A_Rvalue = sqrt(1-exp(-m.A_Deviance/m.T));
    B_Rvalue = sqrt(1-exp(-m.B_Deviance/m.T));

    % plot the coefficient of variation as a matrix (only if significant)
    subplot(3,2,1)
    imagesc(A_Rvalue.*(m.A_pval<o.threshold)); ylabel('Effect on y(t)'); xlabel('Cause y(t-1)');
    title('A R-values'); axis equal; axis tight; xlabel(colorbar,'R');
    clim([0 max(max(A_Rvalue-diag(diag(A_Rvalue))))]); % clip the diagonal values which are much larger

    subplot(3,2,2)
    imagesc(B_Rvalue.*(m.B_pval<o.threshold)); ylabel('Effect on y(t)'); 
    xlabel('Cause x(t)'); set(gca,'xticklabel',{})
    title('B R-values');  axis tight; xlabel(colorbar,'R');

    % show filters as a matrix
    subplot(3,1,2)
    A = m.A; % for i=1:Dy, A(:,i,i)=0; end
    imagesc((1:na)/o.fs,1:Dy,reshape(permute(A,[2 3 1]),Dy,Dy*na)); 
    clim([-1 1]*max(abs(A(:))));   colorbar
    ylabel('Effect on y(t)'); set(gca,'xticklabel',{})
    title('AR filter A')
    [H,Ainv] = varx_trf(m.B,m.A,o.duration*o.fs);
    if 0
        subplot(6,1,4)
        imagesc((1:o.duration*o.fs)/o.fs,1:Dy,reshape(permute(Ainv,[2 3 1]),Dy,Dy*o.duration*o.fs)); 
        clim([-1 1]*max(abs(Ainv(:)))); colorbar
        ylabel('response in y(t)'); set(gca,'xticklabel',{})
        title('Impulse Response of 1/A')
    end
    xlabel('...from y(t-1), lag');

    if o.duration, subplot(3,2,5), else subplot(3,1,3), end
    imagesc((1:size(m.B,1))/o.fs,1:Dy*Dx,reshape(m.B,nb,Dy*Dx)'); ylabel('Effect on y(t)'); 
    clim([-1 1]*max(abs(m.B(:)))); colorbar
    hold on; for i=1:Dx-1, plot([1 size(m.B,1)]/o.fs,[Dy Dy]*i+0.5,'k'); end; hold off
    if o.fs==1, xlabel('... from x(t), lag (samples)'); else xlabel('... from x(t), lag (seconds)'); end
    title('MA filter B')
    if o.duration
        subplot(3,2,6)
        imagesc((1:size(H,1))/o.fs,1:Dy*Dx,reshape(H,o.duration*o.fs,Dy*Dx)'); 
        clim([-1 1]*max(abs(H(:)))); colorbar
        ylabel('response in y(t)'); 
        hold on; for i=1:Dx-1, plot([1 size(H,1)]/o.fs,[Dy Dy]*i+0.5,'k'); end; hold off
        if o.fs==1, xlabel('... from x(t), lag (samples)'); else xlabel('... from x(t), lag (seconds)'); end
        title('Impulse Response of B/A')
    end

end

end

% got this from Bert de Vries's stkplt
function stack

% find plots, get rid of vertical lines
hp = findobj(gca,'type','line'); hplots=[];
for i=1:length(hp)
    xd = get(hp(i),'xdata');
    if any(xd~=xd(1)), hplots=[hplots,hp(i)]; end %if
end
nsig = length(hplots);

grp = 1:nsig; ngrp = max(grp);

space = 0.1/ngrp;  
range = (1-(ngrp+1)*space)/ngrp;

% stack
bias = zeros(nsig,1); scale = zeros(nsig,1);
for i=1:nsig
    y = get(hplots(i),'YData');
    ymax = nanmax(y); ymin = nanmin(y);
    offset = grp(i)*space + (grp(i)-1)*range;
    scale(i) = range/(ymax-ymin + eps); % scale to [0->range]
    bias(i) = offset - scale(i)*ymin;
    set(hplots(i),'YData',scale(i)*y+bias(i), 'Userdata',[scale(i) bias(i)]);
end
set(gca,'YTick',[]);

end
