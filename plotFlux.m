function [] = plotFlux(res, R)
%% plots the flux rates of the specified reaction of the specified simulation result
%
% Input:
%   res:   result structure
%    R:   index of reaction or reaction name
%

if ~exist('res', 'var'), error('Parameter ''res'' is missing.'); end
if ~exist('R', 'var'), disp('Parameter ''R'' is missing. Plotting all...');
    R = 1:length(res.Prob.model.rxns);
end

if ischar(R),
    % reaction name was supplied, transform it into an index
    tmp = find(ismember(res.Prob.model.rxns,R));
    if isempty(tmp), error('Reaction ''%s'' is not known.', R); end
    R = tmp;
end

for iR = R
    
    figure(iR)
    
    % collect x and y data for the plot
    % xdata is a column vector with the middle of the intervals
    % ydata is a 2D-matrix with 3 columns, the first column belongs to the
    % flux rates in x0 solution, the 2nd and the 3rd column belong to
    % the lower and upper bounds. The rows correspond to rows in
    % xdata
    
    xdata = cumsum(res.Prob.model.dT)-res.Prob.model.dT/2;
    ydata = [arrayfun(@(iT) res.Prob.x0(res.Prob.Vars.Fluxes(iT).Indices(iR)), 1:res.Prob.model.nT)'  ...
        reshape(res.FluxVar.Values(iR,:,:),res.Prob.model.nT,2)];
    
    % plot anaerobic/aerobic areas
    switch_time_idx = find(res.Prob.model.ub_var.ETC,1);
    switch_time = xdata(switch_time_idx);
    half_interval_time_correction = res.Prob.model.T/res.Prob.model.nT/2;
    yl=[-100 100];
    area([xdata(1) switch_time]-half_interval_time_correction,yl(2)*ones(1,2),yl(1),'FaceAlpha',0.15,'FaceColor','y','EdgeColor','none'); hold on
    area([switch_time xdata(end)]-half_interval_time_correction,yl(2)*ones(1,2),yl(1),'FaceAlpha',0.1,'FaceColor','m','EdgeColor','none'); hold on
    
    % plot lower bound to upper bound area
    fill([xdata; xdata(end:-1:1)],[ydata(:,end); ydata(end:-1:1,2)],[1 1 1]*0.7, 'EdgeColor','none');
    set(gca, 'FontSize', 14, 'FontName', 'Arial Narrow');
    hold on;
    % plot x0 solution
    plot(xdata,ydata(:,1),'k-o','LineWidth',2,'MarkerFaceColor','k','MarkerSize',3)%,'Clipping','off');
    hold off;
    YLim = get(gca, 'YLim');
    % increase of y-axes length if curves are too close to upper y-axes limit
    %YLim(2) = max(YLim(2), 1.2*max(ydata(:))-YLim(1)*0.2);
    YLim(1) = -0.4;
    YLim(2) = 0.6;
    
    set(gca, 'YLim', YLim);
    set(gca, 'XLim', [0 res.Prob.model.T-0.25]);
    if res.Prob.model.T==24, set(gca,'XTick', 0:6:24); end
    set(gca,'XTick', 0:1:5);
    ylabel({strcat(strrep(res.Prob.model.rxns{iR},'_','\_'),' rate'), '[mmol/g_{DW0}/h]'});
    xlabel('Time [h]');
    % Make axes go through origin instead of left and bottom sides of axes box.
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    
    %NAME = [res.Prob.model.rxns{iR} '_' res.Prob.model.ID]
    
    %saveas(gcf, [NAME '.fig' ] , 'fig')
    %saveas(gcf, [NAME '.eps' ]  , 'epsc')
    
end
end