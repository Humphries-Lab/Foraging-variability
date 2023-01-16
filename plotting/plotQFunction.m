%% plotQFunction - THIS NEEDS WORK

fig = figure;
if max(gT_Q_leave_model.state) > 1 % if model with states for Q leave, include 2 more subplots for Q leave
    tiledlayout(2,2)
    fig.Position = [100 100 700 700];
else
    tiledlayout(1,2)
end

gT_Q_stay_model= gT_Q_stay_model(gT_Q_stay_model.state <=15, :); % only take first states 

if max(gT_Q_leave_model.state) > 6  % for off policy models where extra state at end of leave, remove for plotting purposes
    gT_Q_leave_model= gT_Q_leave_model(gT_Q_leave_model.state <=6,:);
end

tmp_stay_poor_l = table2array(gT_Q_stay_model(gT_Q_stay_model.environ_type==2 & gT_Q_stay_model.patch_type==1, "mean_Q_stay"));
tmp_stay_poor_m = table2array(gT_Q_stay_model(gT_Q_stay_model.environ_type==2 & gT_Q_stay_model.patch_type==2, "mean_Q_stay"));
tmp_stay_poor_h = table2array(gT_Q_stay_model(gT_Q_stay_model.environ_type==2 & gT_Q_stay_model.patch_type==3, "mean_Q_stay"));

% plot Q_stay - poor
ax1 = nexttile;
plot(1:size(tmp_stay_poor_l,1), tmp_stay_poor_l, 1:size(tmp_stay_poor_m,1), tmp_stay_poor_m, 1:size(tmp_stay_poor_h,1), tmp_stay_poor_h);
xlim([1,max([size(tmp_stay_poor_l, 1), size(tmp_stay_poor_m, 1), size(tmp_stay_poor_h, 1)])]);
ax = gca;
ax.FontSize = 14;
xlabel('State (timestep)','FontSize', 16, 'FontName', 'Helvetica');
ylabel('Q-stay','FontSize', 16, 'FontName', 'Helvetica');
set(findobj(gca,'type','line'),'linew',2)
subtitle('Poor environment')

tmp_stay_rich_l = table2array(gT_Q_stay_model(gT_Q_stay_model.environ_type==1 & gT_Q_stay_model.patch_type==1, "mean_Q_stay"));
tmp_stay_rich_m = table2array(gT_Q_stay_model(gT_Q_stay_model.environ_type==1 & gT_Q_stay_model.patch_type==2, "mean_Q_stay"));
tmp_stay_rich_h = table2array(gT_Q_stay_model(gT_Q_stay_model.environ_type==1 & gT_Q_stay_model.patch_type==3, "mean_Q_stay"));

% plot Q_stay - rich
ax2 = nexttile;
plot(1:size(tmp_stay_rich_l,1), tmp_stay_rich_l, 1:size(tmp_stay_rich_m,1), tmp_stay_rich_m, 1:size(tmp_stay_rich_h,1), tmp_stay_rich_h);
xlim([1,max([size(tmp_stay_rich_l, 1), size(tmp_stay_rich_m, 1), size(tmp_stay_rich_h, 1)])]);
ax = gca;
ax.FontSize = 14;
xlabel('State (timestep)','FontSize', 16, 'FontName', 'Helvetica');
ylabel('Q-stay','FontSize', 16, 'FontName', 'Helvetica');
set(findobj(gca,'type','line'),'linew',2)
subtitle('Rich environment')

linkaxes([ax1, ax2], 'y')

% plot Q_leave if multiple states
if max(gT_Q_leave_model.state) > 1

    tmp_leave_poor = [table2array(gT_Q_leave_model(gT_Q_leave_model.environ_type==2 & gT_Q_leave_model.patch_type==1, "mean_Q_leave")),...
        table2array(gT_Q_leave_model(gT_Q_leave_model.environ_type==2 & gT_Q_leave_model.patch_type==2, "mean_Q_leave")),...
        table2array(gT_Q_leave_model(gT_Q_leave_model.environ_type==2 & gT_Q_leave_model.patch_type==3, "mean_Q_leave"))];

    ax3 = nexttile;
    plot(1:max(gT_Q_leave_model.state), tmp_leave_poor);
    xlim([1, max(gT_Q_leave_model.state)])
    ax = gca;
    ax.FontSize = 14;
    xlabel('State (timestep)','FontSize', 16, 'FontName', 'Helvetica');
    ylabel('Q-leave','FontSize', 16, 'FontName', 'Helvetica');
    set(findobj(gca,'type','line'),'linew',2)
    subtitle('Poor environment')

    tmp_leave_rich = [table2array(gT_Q_leave_model(gT_Q_leave_model.environ_type==1 & gT_Q_leave_model.patch_type==1, "mean_Q_leave")),...
        table2array(gT_Q_leave_model(gT_Q_leave_model.environ_type==1 & gT_Q_leave_model.patch_type==2, "mean_Q_leave")),...
        table2array(gT_Q_leave_model(gT_Q_leave_model.environ_type==1 & gT_Q_leave_model.patch_type==3, "mean_Q_leave"))];

    ax4 = nexttile;
    plot(1:max(gT_Q_leave_model.state), tmp_leave_rich);
    xlim([1, max(gT_Q_leave_model.state)])
    ax = gca;
    ax.FontSize = 14;
    xlabel('State (timestep)','FontSize', 16, 'FontName', 'Helvetica');
    ylabel('Q-leave','FontSize', 16, 'FontName', 'Helvetica');
    set(findobj(gca,'type','line'),'linew',2)
    subtitle('Rich environment')

    linkaxes([ax3, ax4], 'y')
end

leg = legend('Low', 'Medium', 'High');
leg.Layout.Tile = 'north';
title(leg, 'Patch type')
sgtitle(sprintf('State-action values by patch and environment. M%d.', model),'FontSize', 16, 'FontName', 'Helvetica', 'FontWeight', 'Bold');
