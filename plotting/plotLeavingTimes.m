function plotLeavingTimes(model)
%% figure - plotLeavingTimes

load('~/Dropbox/foraging/raw_data/summary/young_variables/subMean_young.mat') % load real data

%%Sub-plot
% model = 1:26; % if want to plot more than one model, uncomment

% set up colours
c = [0.6353    0.0784    0.1843; 0    0.4471    0.7412]; % [rich, poor]

%% plot the actual experimental data
NSub = size(subMean_young, 3);
SubLT(:,:,1) = squeeze(subMean_young(1,:,:)); % rich environment
SubLT(:,:,2) = squeeze(subMean_young(2,:,:)); % poor environment
SubLTMean = [mean(SubLT(:,:,1), 2), mean(SubLT(:,:,2), 2)];
SubSEM = [std(SubLT(:,:,1),[],2)./sqrt(NSub), std(SubLT(:,:,2), [], 2)./sqrt(NSub)]; % Standard Error

ts = tinv([0.025  0.975],NSub-1);      % T-Score
SubCI = ts(2).*SubSEM;

figure()
t = tiledlayout('flow', 'TileSpacing', 'Compact');
title(t,'Mean patch leaving times');
xlabel(t, 'Patch yield type')
ylabel(t, 'Time in patch (s)')

for m = model
    nexttile
    errorbar([1 2 3],SubLTMean(:,1), SubCI(:,1),'LineStyle','--', 'Color',c(1,:) ,'LineWidth', 1);
    hold on
    errorbar([1 2 3],SubLTMean(:,2), SubCI(:,2),'LineStyle','--', 'Color',c(2,:) ,'LineWidth', 1);

    %% overlay simulated data
    load(sprintf('~/Dropbox/foraging/outputs/M%d/SimulationResults.mat', model)) % load simulated results for each model

    ts = tinv([0.025  0.975],NSim-1);      % T-Score
    LTCI = LTSEM .* ts(2);

    errorbar([1 2 3], LTMean(:,1), LTCI(:,1), 'LineWidth', 2, 'Color', c(1,:)); % rich environment
    errorbar([1 2 3], LTMean(:,2), LTCI(:,2), 'LineWidth', 2, 'Color', c(2,:)); % poor environment
    hold off

    title(sprintf('M%d', model))
    % tidy up axes
    ax = gca;
    axis(ax, 'tight')
    xlim(ax, xlim(ax) + [-1,1]*range(xlim(ax)).* 0.1)
    xticks(ax, [1 2 3])
    xticklabels(ax, {'Low'; 'Medium'; 'High'})
    ylim([0, 30])

end

leg = legend('Rich environment', 'Poor environment');
leg.Layout.Tile = 'north';
subtitle(t, sprintf('Error bars show 95%% CI. NSim = %d', NSim))

