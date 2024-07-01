% got this from Bert de Vries's stkplt
function stackplot

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