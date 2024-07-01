function [Graph,G_plot]=varx_display(m,o)
% [Graph,G_plot]=varx_display(model,options)
% Displays the model effect sizes and model filter parameters. Example:
%
% model = varx(y,na,x,nb);
% varx_display(model,xname={},yname={},duration=0,fs=1,threshold=0.001,plottype='Default');
%
% The argument model is a structure with the following stucture elements
% A,B,A_pval,B_pval,A_Deviance,B_Deviance,T.base as created by model =
% varx().
%
% options are name=value optional arguments: yname,xname,duration,fs,threshold, plottype
%
% yname,xname: cell arrays that specify the names of input x and output
%    y in the model. If they are provided, and plottype 'Graph' is
%    selected, then model structure is displayed as a circular graph, and
%    Graph and G_plot can be used to plot with:
%
% plot(Graph,'LineWidth',G_plot.width,'XData',G_plot.xdata,'YData',G_plot.ydata,...
%            'EdgeColor',G_plot.color,'NodeColor',G_plot.nodecolor);
%
% plottype: if yname is not provided, plottype defaults to 'Default'
%    'Default': everything is shown as matrix. Suggest to use this
%           when there are a lot of variables in the model.
%    'Matrix':
%           matrices with input x and output y names.
%    'Graph': digraph visualization, Effect size matrices with input x
%           and output y names, and impulse responses, also
%           outputs digraph structures. You can add 'no filter' to the
%           string if you want to suppress showing the filters. 
%
% duration: a scalar indicating for how long to compute and display the
%   impulse respones (IR). The IR are from every variable to every
%   other variable. If it is omited then the IR is not displayed.
%
% fs: the sampling rate. If given, then filters and impulse response are
%   displayed in seconds.
%
% threshold: Only links with pvalue below this threshold will be displayed
%    (defaults to 0.001)

% March 11, 2024, Lucas Parra
% March 16, 2024, changed to display R as effect size. added threshold as argument, and added options arguments
% April 1, 2024, Aimar Silvan, added display plottype: 'Matrix', 'Graph',
%                              'Default', now also plots R-value matrices
% April 11, 2024, Lucas, fix impulse response display in 'Default' option
% April 13, 2024, Lucas, make more space for Graph, and change a few other visuals

% define arguments with defaults for the options variables o
arguments
    m
    o.yname cell = {}     % show as matrix
    o.xname cell = {}
    o.duration double = 0 % no impulse response is shown
    o.fs double = 1       % lags displayed in samples.
    o.threshold double = 0.001 % only links with lower values shown
    o.plottype char = 'Default'
end

clf

[nb,~,Dx] = size(m.B);
[na,~,Dy] = size(m.A);

if ~isempty(o.yname) && contains(o.plottype, 'Graph')


    % combine pval and effect size for single graph for A and B
    pval   = [m.A_pval   m.B_pval;   ones(Dx,Dy+Dx)];
    Rvalue = [m.A_Rvalue m.B_Rvalue; zeros(Dx,Dy+Dx)];

    %tiledlayout(4,Dx+Dy)%    subplot(2,3,1)
    %nexttile([2 floor((Dx+Dy)/2)])
    subplot(2,2,1)

    % plot connectivity as a graph
    node_name = [o.yname(:)' o.xname(1:Dx)];

    % set graph visuals
    h=plot(digraph(zeros(Dy)),'Layout','circle');
    x = linspace(.8,-0.8,Dx); x=x-mean(x);
    G_plot.ydata=[get(h,'YData') x];
    G_plot.xdata=[get(h,'XData') -2*ones(1,Dx)];
    G_plot.nodecolor = [repmat([0 0 0.75],Dy,1); repmat([0.75 0 0],Dx,1)];

    % generate the graph and set more graph visuals
    [s,t]=meshgrid(1:Dy+Dx,1:Dy+Dx);
    Adj = pval<o.threshold;
    Graph = digraph(s(Adj),t(Adj),Rvalue(Adj),node_name,'omitselfloops');
    G_plot.width = 10*Graph.Edges.Weight/max(Graph.Edges.Weight);
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
    plot(Graph,'LineWidth',G_plot.width,'XData',G_plot.xdata,'YData',G_plot.ydata,'EdgeColor',G_plot.color,'Nodecolor',G_plot.nodecolor,'ArrowSize',20);
    ylim([-1 1]); axis equal; axis off; 
    A_Rvalue = m.A_Rvalue;
    B_Rvalue = m.B_Rvalue;

    % plot the coefficient of variation as a matrix (only if significant)
    A_Adj = Adj(1:length(o.yname),1:length(o.yname));
    subplot(2,4,3) % nexttile([2 floor((Dx+Dy)/3)])
    imagesc(A_Rvalue.*A_Adj); ylabel('Effect on y(t)'); xlabel('Cause y(t-1)');
    clim_max = max(max(A_Rvalue-diag(diag(A_Rvalue)))); % clip the diagonal values which are much larger
    clim_max = max([clim_max; B_Rvalue(:)]); % use the same clim for A and B
    if clim_max>0, clim([0 clim_max]); end
    title('A Effect size');
    set(gca,"YTick",[1:length(o.yname)]); set(gca,"YTickLabel",node_name(1:length(o.yname)),"YTickLabelRotation",45)
    set(gca,"XTick",[1:length(o.yname)]); set(gca,"XTickLabel",node_name(1:length(o.yname)),"XTickLabelRotation",45)


    if ~isempty(B_Rvalue)
        % nexttile([2 Dx+Dy-floor((Dx+Dy)/2)-floor((Dx+Dy)/3)])
        subplot(2,4,4)
        imagesc(B_Rvalue.*(m.B_pval<o.threshold)); ylabel('Effect on y(t)');
        if clim_max>0, clim([0 clim_max]); end
        xlabel('Cause x(t)'); set(gca,'xticklabel',{})
        title('B Effect size');
        set(gca,'DataAspectRatio',[1 1 1])
        %      set(gca,"YTick",[1:size(A_Rvalue,1)]); set(gca,"YTickLabel",node_name(1:size(A_Rvalue,2)),"YTickLabelRotation",45)
        set(gca,"XTick",[1:size(B_Rvalue,1)]); set(gca,"XTickLabel",node_name(end-size(B_Rvalue,2)+1:end),"XTickLabelRotation",45)
    end
    xlabel(colorbar,'R');


    if ~contains(o.plottype, 'no filter') % unless the user wants no filters ...
        % make and plot the filters
        if na>1 % it is na is less, showing as a curve is not helpful
            for i=1:Dy
                subplot(4,Dy+Dx,2*(Dy+Dx)+i);
                plot((1:na)/o.fs,m.A(:,:,i).*shiftdim(pval(1:Dy,i)<0.0001,1)); stackplot
                set(gca,'ytick',[]); title(['A: ' node_name{i}])
                if o.fs==1, xlabel('samples'); else, xlabel('seconds'); end
                if i==1
                    axis on; ylabel('Filters')
                    pos=get(gca,'Position');
                    h=legend(node_name{1:Dy},'Location','westoutside');
                    set(gca,'Position',pos);
                    title(h,'Effect on')
                end
            end
            for i=1:Dx
                subplot(4,Dy+Dx,2*(Dy+Dx)+Dy+i);
                plot((1:nb)/o.fs*1,m.B(:,:,i).*shiftdim(pval(1:Dy,Dy+i)<0.0001,1));
                set(gca,'ytick',[]); title(['B: ' node_name{Dy+i}]); stackplot
                if o.fs==1, xlabel('samples'); else, xlabel('seconds'); end
            end
        else  % if na<=1 better show filters as a matrix
            if Dx>0, subplot(4,2,5); else subplot(4,1,3); end;
            A = m.A; % for i=1:Dy, A(:,i,i)=0; end
            imagesc((1:na*Dy)/o.fs,1:Dy,reshape(permute(A,[2 3 1]),Dy,Dy*na));
            clim([-1 1]*max(abs(A(:))));   colorbar
            ylabel('Effect on y(t)');
            title('Filter A')

            if Dx>0
                subplot(4,2,6);
                imagesc((1:nb)/o.fs,1:Dy*Dx,reshape(m.B,nb,Dy*Dx)'); ylabel('Effect on y(t)');
                if ~isempty(max(abs(m.B(:))))
                    clim([-1 1]*max(abs(m.B(:))));
                end
                hold on; for i=1:Dx-1, plot([1 size(m.B,1)]/o.fs,[Dy Dy]*i+0.5,'k'); end; hold off
                title('Filter B')
            end

        end

    end

    % show the impulse responses for duration
    if o.duration
        [H,Ainv] = varx_trf(m.B,m.A,o.duration*o.fs);
        for i=1:Dy
            subplot(4,Dy+Dx,3*(Dy+Dx)+i);
            plot((1:o.duration*o.fs)/o.fs,Ainv(:,:,i));
            set(gca,'ytick',[]); stackplot
            if o.fs==1, xlabel('samples'); else, xlabel('seconds'); end
            if i==1
                axis on; ylabel('Impulse response')
                pos=get(gca,'Position');
                h=legend(node_name{1:Dy},'Location','westoutside');
                set(gca,'Position',pos);
                title(h,'Response of')
            end
        end
        for i=1:Dx
            subplot(4,Dy+Dx,3*(Dy+Dx)+Dy+i);
            plot((1:o.duration*o.fs)/o.fs,H(:,:,i));
            set(gca,'ytick',[]); stackplot
            if o.fs==1, xlabel('samples'); else, xlabel('seconds'); end
        end
    end

elseif ~isempty(o.yname) && strcmpi(o.plottype, 'Matrix')

    node_name = [o.yname(:)' o.xname(1:Dx)];
    % sqrt of Coeficient of variation, i.e. R-square
    A_Rvalue = m.A_Rvalue;
    B_Rvalue = m.B_Rvalue;

    % plot the coefficient of variation as a matrix (only if significant)
    t = tiledlayout(1,2);
    nexttile
    imagesc(A_Rvalue.*(m.A_pval<o.threshold)); ylabel('Effect on y(t)'); xlabel('Cause y(t-1)');
    title('A Effect size'); axis equal; axis tight;
    clim_max = max(max(A_Rvalue-diag(diag(A_Rvalue)))); % clip the diagonal values which are much larger
    clim_max = max([clim_max; B_Rvalue(:)]); % use the same clim for A and B
    if clim_max>0, clim([0 clim_max]); end
    set(gca,"YTick",[1:size(A_Rvalue,1)]); set(gca,"YTickLabel",node_name(1:size(A_Rvalue,2)),"YTickLabelRotation",45)
    set(gca,"XTick",[1:size(A_Rvalue,1)]); set(gca,"XTickLabel",node_name(1:size(A_Rvalue,2)),"XTickLabelRotation",45)
    if ~isempty(B_Rvalue)
        nexttile
        imagesc(B_Rvalue.*(m.B_pval<o.threshold)); ylabel('Effect on y(t)');
        if clim_max>0, clim([0 clim_max]); end
        xlabel('Cause x(t)'); axis equal; axis tight;set(gca,'xticklabel',{})
        title('B Effect size');  axis tight; xlabel(colorbar,'R');
        set(gca,"YTick",[1:size(A_Rvalue,1)]); set(gca,"YTickLabel",node_name(1:size(A_Rvalue,2)),"YTickLabelRotation",45)
        set(gca,"XTick",[1:size(B_Rvalue,1)]); set(gca,"XTickLabel",node_name(end-size(B_Rvalue,2)+1:end),"XTickLabelRotation",45)
    end
    xlabel(colorbar,'R');

elseif isempty(o.yname) || strcmpi(o.plottype, 'Default')

    if o.duration, rows=4; else rows=3; end

    % sqrt of Coeficient of variation, i.e. R-square
    A_Rvalue = m.A_Rvalue;
    B_Rvalue = m.B_Rvalue;

    % plot the coefficient of variation as a matrix (only if significant)
    tiledlayout(rows,3)
    nexttile([1 2])
    imagesc(A_Rvalue.*(m.A_pval<o.threshold)); ylabel('Effect on y(t)'); xlabel('Cause y(t-1)');
    title('A Effect size'); axis equal; axis tight; xlabel(colorbar,'R');
    clim([0 max(max(A_Rvalue-diag(diag(A_Rvalue))))]); % clip the diagonal values which are much larger

    nexttile
    imagesc(B_Rvalue.*(m.B_pval<o.threshold)); ylabel('Effect on y(t)');
    xlabel('Cause x(t)'); axis equal; axis tight;set(gca,'xticklabel',{})
    title('B Effect size');  axis tight; xlabel(colorbar,'R');

    % show filters as a matrix
    nexttile([1 3])
    A = m.A; % for i=1:Dy, A(:,i,i)=0; end
    imagesc(0.5:na+0.5,1:Dy,reshape(permute(A,[2 3 1]),Dy,Dy*na));
    set(gca,'XTick',1:na)
    clim([-1 1]*max(abs(A(:))));   colorbar
    ylabel('Effect on y(t)');
    xlabel('lag (samples)');
    title('Filter A')

    nexttile([1 3])
    imagesc((1:nb)/o.fs,1:Dy*Dx,reshape(m.B,nb,Dy*Dx)'); ylabel('Effect on y(t)');
    if ~isempty(max(abs(m.B(:))))
        clim([-1 1]*max(abs(m.B(:))));
    end
    colorbar
    hold on; for i=1:Dx-1, plot([1 size(m.B,1)]/o.fs,[Dy Dy]*i+0.5,'k'); end; hold off
    title('Filter B')

    if o.duration
        nexttile([1 3])
        [H,Ainv] = varx_trf(m.B,m.A,o.duration*o.fs);
        imagesc((1:size(H,1))/o.fs,1:Dy*Dx,reshape(H,o.duration*o.fs,Dy*Dx)');
        clim([-1 1]*max(abs(H(:)))); colorbar
        ylabel('response in y(t)');
        hold on; for i=1:Dx-1, plot([1 size(H,1)]/o.fs,[Dy Dy]*i+0.5,'k'); end; hold off
        title('Impulse Response H=B/(1-A)')
    end
    if o.fs==1, xlabel('lag (samples)'); else xlabel('lag (seconds)'); end

end

end


