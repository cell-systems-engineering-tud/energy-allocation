function [] = plotCompound(res, M)
%% plots the amounts of the specified compound of the specified simulation result
%
% Input:
%   res:   result structure
%    M:   index of imbalanced compound or compound name
%

if ~exist('res', 'var'), error('Parameter ''res'' is missing.'); end
if ~exist('M', 'var'), disp('Parameter ''iM'' is missing. Plotting all...');
M = 1:length(res.Prob.model.ImbalancedMets);
end

if ischar(M),
    % compound name was supplied, transform it into an index
    tmp = find(ismember(res.Prob.model.mets(res.Prob.model.ImbalancedMets),M));
    if isempty(tmp), error('Compound ''%s'' is not known.', M); end
    M = tmp;
end


for iM = M
    figure(iM)
    
    % collect x and y data for the plot
    % xdata is a column vector with the time points
    % ydata is a 2D-matrix with 3 columns, the first column belongs to the
    % compound amounts in x0 solution, the 2nd and the 3rd column belong to
    % the lower and upper bounds. The rows correspond to the time points in
    % xdata
    
    xdata = [0; cumsum(res.Prob.model.dT)];
    ydata = [arrayfun(@(iT) res.Prob.Vars.ImbMets.S(iM,1:res.Prob.Vars.ImbMets.Tind(iT))*res.Prob.x0(1:res.Prob.Vars.ImbMets.Tind(iT)), 1:res.Prob.model.nT)'  ...
        reshape(res.CompoundVar.Values(iM,:,:),res.Prob.model.nT,2)];
    ydata = [[res.Prob.x0(res.Prob.Vars.ImbMets.StartAmounts(iM)) res.CompoundStartVar.Values(iM,:)]; ydata];
    
    % plot anaerobic/aerobic areas
    switch_time_idx = find(res.Prob.model.ub_var.ETC,1);
    switch_time = xdata(switch_time_idx);
    yl=[0 10];
    area([xdata(1) switch_time],yl(2)*ones(1,2),'FaceAlpha',0.15,'FaceColor','y','EdgeColor','none'); hold on
    area([switch_time xdata(end)],yl(2)*ones(1,2),'FaceAlpha',0.1,'FaceColor','m','EdgeColor','none'); hold on
    
    % plot lower bound to upper bound area
    fill([xdata; xdata(end:-1:1); xdata(1)],[ydata(:,end); ydata(end:-1:1,2); ydata(1,end)],[1 1 1]*0.7, 'EdgeColor','none');
    %   area(xdata, ydata(:,end),'FaceColor',[1 1 1]*0.7, 'EdgeColor','none'); hold on
    %   area(xdata, ydata(:,2),'FaceColor','w', 'EdgeColor','none');
    
    set(gca, 'FontSize', 14, 'FontName', 'Arial Narrow');
    hold on;
    % plot x0 solution
    plot(xdata,ydata(:,1),'k-o','LineWidth',2,'MarkerFaceColor','k','MarkerSize',3,'Clipping','off');
    hold off;
    YLim = get(gca, 'YLim');
    % increase of y-axes length if curves are too close to upper y-axes limit
    YLim = [0 max(max(ydata(:))*1.2,YLim(2)*0.1)];
    %  YLim = [0 max(ydata(:))*1.2];
    
    if YLim(1)>YLim(2)
        YLim = wrev(YLim);
    end
    set(gca, 'YLim', YLim);
    set(gca, 'XLim', [0 res.Prob.model.T]);
    if res.Prob.model.T==24, set(gca,'XTick', 0:6:24); end
    ylabel(sprintf('%s amount [mmol]', strrep(res.Prob.model.mets{res.Prob.model.ImbalancedMets(iM)},'_','\_')));
    xlabel('Time [h]');
end
end