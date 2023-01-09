function plot_leaving_stats_line(gT_model, subMean_young)

close all
% set up colours
c = [0.6353    0.0784    0.1843; 0    0.4471    0.7412]; % [rich, poor]

%% plot the actual experimental data
subMean_poor = squeeze(subMean_young(2,:,:))';
subMean_rich = squeeze(subMean_young(1,:,:))';
mean_poor_LTs = mean(subMean_poor);
mean_rich_LTs = mean(subMean_rich);
SEM_poor = std(subMean_poor)/sqrt(length(subMean_poor));               % Standard Error
SEM_rich = std(subMean_rich)/sqrt(length(subMean_rich));               % Standard Error
ts = tinv([0.025  0.975],length(subMean_poor)-1);      % T-Score
for i = 1:3
    error_poor(i) = 1.96.*SEM_poor(i);
    error_rich(i) = 1.96.*SEM_rich(i);
end

figure()
errorbar([1 2 3],mean_poor_LTs,error_poor,'LineStyle','--', 'Color',c(2,:) ,'LineWidth', 1);
hold on
errorbar([1 2 3], mean_rich_LTs, error_rich,'LineStyle', '--', 'Color',c(1,:),'LineWidth',1)

%% overlay mean leaving times by patch and environment
errorbar(unique(gT_model.patch_type), gT_model.g_mean_lt(4:6), gT_model.g_mean_lt(4:6)-gT_model.g_ci_lt(4:6,1), 'LineWidth', 2, 'Color', c(2,:)); % poor environment
errorbar(unique(gT_model.patch_type), gT_model.g_mean_lt(1:3), gT_model.g_mean_lt(1:3)-gT_model.g_ci_lt(1:3,1), 'LineWidth', 2, 'Color', c(1,:)); % rich environment

hold off

ylim([0 30]);
xlim([0.5 3.5]);
ax=gca;
ax.XTick= [1 2 3];ax.XTickLabel={'Low','Medium','High'};

xlabel('Patch yield type')
ylabel('Time in patch (s)')
legend('Poor environment - Observed', 'Rich environment - Observed', 'Poor environment - Simulated', 'Rich environment - Simulated','Location','northwest');
title(sprintf('M%d', gT_model.model(1)));
subtitle(sprintf('Error bars show 95%% CI. N = %d', max(gT_model.GroupCount)))

%% plot mean reward rates at leaving by patch and environment
% figure('Name', 'Leaving reward rate');
% % errorbar(unique(gT_model.patch_type), gT_model.g_mean_lrr(4:6), gT_model.g_mean_lrr(4:6)-gT_model.g_ci_lrr(4:6,1), 'LineWidth', 2, 'Color', c(2,:));
% errorbar(unique(gT_model.patch_type), gT_model.g_mean_lrr(1:3), gT_model.g_mean_lrr(1:3)-gT_model.g_ci_lrr(1:3,1), 'LineWidth', 2, 'Color', c(1,:));
%
% % tidy up axes
% ax = gca;
% axis(ax, 'tight')
% xlim(ax, xlim(ax) + [-1,1]*range(xlim(ax)).* 0.1)
% xticks(ax, [1 2 3])
% xticklabels(ax, {'Low'; 'Medium'; 'High'})
% ylim([0, 50])
%
% title(sprintf('M%d', gT_model.model(1)));
% xlabel('Patch yield type')
% ylabel('Foreground reward rate')
% legend('Poor environment', 'Rich environment','Location','northwest');
% subtitle(sprintf('Error bars show 95%% CI. N = %d', max(gT_model.GroupCount)))
%
% if save_plots == 1
%     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])
%    print(sprintf('/Users/exs165/Dropbox/foraging-project/results/M%d/M%d_leaving_rr', gT_model.model(1), gT_model.model(1)), '-r600','-dpng')
%
% end
%
% end