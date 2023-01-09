function subplot_leaving_stats_line(gT)

close all
% set up colours
c = [0.6353    0.0784    0.1843; 0    0.4471    0.7412];

%% plot mean leaving times by patch and environment
figure('Name', 'Leaving times');
t = tiledlayout('flow', 'TileSpacing', 'Compact');
title(t,'Mean time in patch by environment and patch type');
xlabel(t,'Patch yield type')
ylabel(t,'Time in patch (s)')

for model = 1:12
    nexttile;
    errorbar(unique(gT{model}.patch_type), gT{model}.g_mean_lt(4:6), gT{model}.g_mean_lt(4:6)-gT{model}.g_ci_lt(4:6,1), 'LineWidth', 1, 'Color', c(2,:));
    hold on 
    errorbar(unique(gT{model}.patch_type), gT{model}.g_mean_lt(1:3), gT{model}.g_mean_lt(1:3)-gT{model}.g_ci_lt(1:3,1), 'LineWidth', 1, 'Color', c(1,:));
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

leg = legend('Poor environment', 'Rich environment');
leg.Layout.Tile = 'north';
subtitle(t, sprintf('Error bars show 95%% CI. N = %d', max(gT{model}.GroupCount)))

if save_plots == 1
    print('/Users/exs165/Dropbox/foraging-project/results/subplot_leaving_times', '-dpdf','-fillpage')
end


%% plot mean reward rates at leaving by patch and environment
figure('Name', 'Leaving reward rate');
t = tiledlayout('flow', 'TileSpacing', 'Compact');
title(t,'Mean reward rate at leaving by environment and patch type');
xlabel(t,'Patch yield type')
ylabel(t,'Foreground reward rate')

for model = 1:12
    nexttile;

    errorbar(unique(gT{model}.patch_type), gT{model}.g_mean_lrr(4:6), gT{model}.g_mean_lrr(4:6)-gT{model}.g_ci_lrr(4:6,1), 'LineWidth', 1, 'Color', c(2,:));
    hold on
    errorbar(unique(gT{model}.patch_type), gT{model}.g_mean_lrr(1:3), gT{model}.g_mean_lrr(1:3)-gT{model}.g_ci_lrr(1:3,1), 'LineWidth', 1, 'Color', c(1,:));
    hold off

    title(sprintf('M%d', model))
    % tidy up axes
    ax = gca;
    axis(ax, 'tight')
    xlim(ax, xlim(ax) + [-1,1]*range(xlim(ax)).* 0.1)
    xticks(ax, [1 2 3])
    xticklabels(ax, {'Low'; 'Medium'; 'High'})
    ylim([0, 55])

end

leg = legend('Poor environment', 'Rich environment');
leg.Layout.Tile = 'north';
subtitle(t, sprintf('Error bars show 95%% CI. N = %d', max(gT{model}.GroupCount)))

end