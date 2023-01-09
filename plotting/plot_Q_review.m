function plot_Q_review(gT_Q_stay_model, model, n_sim, save_plots)

close all
figure

gT_Q_stay_model= gT_Q_stay_model(gT_Q_stay_model.state <=15, :); % only take first states

tmp_stay_rich_l = table2array(gT_Q_stay_model(gT_Q_stay_model.environ_type==2 & gT_Q_stay_model.patch_type==1, "mean_Q_stay"));
tmp_stay_rich_m = table2array(gT_Q_stay_model(gT_Q_stay_model.environ_type==2 & gT_Q_stay_model.patch_type==2, "mean_Q_stay"));
tmp_stay_rich_h = table2array(gT_Q_stay_model(gT_Q_stay_model.environ_type==2 & gT_Q_stay_model.patch_type==3, "mean_Q_stay"));

% plot Q_stay - rich
plot(1:size(tmp_stay_rich_l,1), tmp_stay_rich_l, 1:size(tmp_stay_rich_m,1), tmp_stay_rich_m, 1:size(tmp_stay_rich_h,1), tmp_stay_rich_h);
xlim([1,max([size(tmp_stay_rich_l, 1), size(tmp_stay_rich_m, 1), size(tmp_stay_rich_h, 1)])]);
ax = gca;
ax.FontSize = 14;
xlabel('State (timestep)','FontSize', 16, 'FontName', 'Helvetica');
ylabel('Q-stay','FontSize', 16, 'FontName', 'Helvetica');
set(findobj(gca,'type','line'),'linew',2)
subtitle(sprintf('M%d', model-14))

leg = legend('Low yield', 'Medium yield', 'High yield', 'Location', 'northeast');
title(leg, 'Patch type')

if save_plots == 1
    print(sprintf('/Users/exs165/Dropbox/9-month-review/Images/M%d_Q_values', model), '-dsvg')
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])
    print(sprintf('/Users/exs165/Dropbox/9-month-review/Images/M%d_Q_values', model), '-r600', '-dpng')
end