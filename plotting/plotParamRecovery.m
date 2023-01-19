function plotParamRecovery(SimParams, FitParams)

%% plot correlation - simulated value versus actual value

ParamNames = {'AlphaPatch-Q', 'AlphaRho', 'InverseTemp'};
figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');

% for each parameter, plot sim vs fit 
for i= 1:size(SimParams,1)
    ax = nexttile;
    plot(SimParams(i,:), FitParams(i,:), 'o', 'markersize', 8, 'linewidth', 1)
    if i == 3
    corrcoef(log(SimParams(i,:)),log(FitParams(i,:))) % do correlation between logs of beta's 
    else
    corrcoef(SimParams(i,:),FitParams(i,:))
    end
    %corr(SimParams(i,:)',FitParams(i,:)', 'Type', 'Spearman')
    title(ParamNames{i})
    xlabel(sprintf('Simulated %s', ParamNames{i}))
    ylabel(sprintf('Fit %s', ParamNames{i}))
end
title(tl, 'Parameter recovery')
set(ax,'xscale', 'log', 'yscale' ,'log')

% plot trade off between parameters

figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');
nexttile;
plot(FitParams(1,:), FitParams(2,:), 'o', 'markersize', 8, 'linewidth', 1)
title('fit alpha Q/patch vs fit alpha rho')
xlabel('fit alpha Q/patch')
ylabel('fit alpha rho')

ax = nexttile;
plot(FitParams(2,:), FitParams(3,:), 'o', 'markersize', 8, 'linewidth', 1)
title('fit alpha rho vs fit beta')
xlabel('fit alpha rho')
ylabel('fit beta')
set(ax,'yscale' ,'log')

ax = nexttile;
plot(FitParams(1,:), FitParams(3,:), 'o', 'markersize', 8, 'linewidth', 1)
title('fit alpha Q/patch vs fit beta')
xlabel('fit alpha Q/patch')
ylabel('fit beta')
set(ax,'yscale' ,'log')
title(tl, 'Parameter trade-offs - correlations between fit parameters')