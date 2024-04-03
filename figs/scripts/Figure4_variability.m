%%% Figure 4 foraging paper 
% Emma Scholey 7 Aug 2023

clearvars; close all

addpath('./plotFunctions')
addpath('../../data/fitting_data_240117/')
addpath('../../model/helperFunctions/')
addpath('../../data/experiment_data/')
addpath('../../data/simulation_data/')
addpath('../../data/analytical_data/')

run figure_properties_foraging.m

study = {'leheron', 'contrerashuerta', 'kane'};

save_figs = 1; % whether to save figures

%% Panel: SD_leave predictions for range of beta 

for s = 1:numel(study)
    task = buildTask(study{s});

    load(sprintf('fitting_results_M27_%s',study{s}), '-mat', 'minNLLFitParams');
    beta_rich = mean(minNLLFitParams.beta_rich);     beta_poor = mean(minNLLFitParams.beta_poor);

    load(sprintf('expectedLT_%s_rangebeta', study{s}), '-mat')

    % Panel: expected SD for range of beta

    SD_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    colororder(color.patch);
    semilogx(data.explore,data.SD_leave,'LineWidth',widths.plot); hold on
    xlabel('Beta (higher = exploit)')
    ylabel('SD_{leave}');
    ylim([0 5])
    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    line([beta_rich beta_rich],[0 SD_fig.Children.YLim(2)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
    line([beta_poor beta_poor],[0 SD_fig.Children.YLim(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')

    if save_figs == 1
        print([export_path, sprintf('fig4_%s_expected_SD', study{s})],'-dsvg')
        print([overleaf_path, 'fig4/', sprintf('fig4_%s_expected_SD', study{s})],'-dsvg')
    end

end


%% Panels: simulated vs empirical mean leaving times
% for winning model, M27
for s = 1:numel(study)
    task = buildTask(study{s});

    load(sprintf('fitting_results_M27_%s',study{s}), '-mat');
    load(sprintf('subjectLT_%s', study{s}), '-mat');
    load(sprintf('modelLT_%s_M27', study{s}), '-mat')

    % summarise data
    modelSD = mean(log(simulated_leave_times.sd),3); modelSEM = std(log(simulated_leave_times.sd), [], 3)./sqrt(size(simulated_leave_times.sd,3));
    subjectSD = mean(log(subject_leave_times.sd),3); subjectSEM = std(log(subject_leave_times.sd), [], 3)./sqrt(size(subject_leave_times.sd,3));

    % SD panel
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.rectangle)
    tl = tiledlayout(1,2, 'TileSpacing', 'tight');
    nexttile
    % plot the actual experimental data
    errorbar(subjectSD(1,:),subjectSEM(1,:),lines.exp, 'Color',color.rich ,'LineWidth', widths.plot); hold on
    errorbar(subjectSD(2,:),subjectSEM(2,:),lines.exp, 'Color',color.poor,'LineWidth', widths.plot)

    set(gca,'XLim',[0, task.nPatch+1],'XTick',1:task.nPatch,'XTickLabel',task.patchNames, 'YLim', [0 6])
    ylabel('SD_{leave} (s)')

    % plot the simulated data on another panel
    nexttile
    errorbar(modelSD(1,:),modelSEM(1,:),lines.model, 'Color',color.rich ,'LineWidth', widths.plot), hold on
    errorbar(modelSD(2,:),modelSEM(2,:),lines.model, 'Color',color.poor,'LineWidth', widths.plot)
    set(gca,'XLim',[0, task.nPatch+1],'XTick',1:task.nPatch,'XTickLabel',task.patchNames, 'YLim', [0 6])
    set(gca, 'YTick', [], 'YColor', 'none');

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig4_%s_simulated_sd_LT', study{s})],'-dsvg')
        print([overleaf_path, 'fig4/', sprintf('fig4_%s_simulated_sd_LT', study{s})],'-dsvg')
    end

    % plot on one panel to show overlap
     % sd panel
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square)
    nexttile
    % plot the actual experimental data
    errorbar(subjectSD(1,:),subjectSEM(1,:),lines.exp, 'Color',color.rich ,'LineWidth', widths.plot); hold on
    errorbar(subjectSD(2,:),subjectSEM(2,:),lines.exp, 'Color',color.poor,'LineWidth', widths.plot)

    set(gca,'XLim',[0, task.nPatch+1],'XTick',1:task.nPatch,'XTickLabel',task.patchNames, 'YLim', [0 6])
    ylabel('SD_{leave} (s)')

    errorbar(modelSD(1,:),modelSEM(1,:),lines.model, 'Color',color.rich ,'LineWidth', widths.plot)
    errorbar(modelSD(2,:),modelSEM(2,:),lines.model, 'Color',color.poor,'LineWidth', widths.plot)

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig4_%s_simulated_sd_LT_overlap', study{s})],'-dsvg')
        print([overleaf_path, 'fig4/', sprintf('fig4_%s_simulated_sd_LT_overlap', study{s})],'-dsvg')
    end

    % write to csv for equivalence testing in R
    tmp = [squeeze(subject_leave_times.sd(1,:,:))', squeeze(subject_leave_times.sd(2,:,:))'];
    writematrix(tmp,sprintf('../../data/experiment_data/figures_stats_data/fig4_%s_subject_SD_data.csv', study{s}))

    tmp = [squeeze(simulated_leave_times.sd(1,:,:))', squeeze(simulated_leave_times.sd(2,:,:))'];
    writematrix(tmp,sprintf('../../data/experiment_data/figures_stats_data/fig4_%s_sim_SD_data.csv', study{s}))

end

%% Panel: per-participant fit (simulated vs actual) 
% for winning model, M27
for s = 1:numel(study)
    task = buildTask(study{s});

    % get sd leave times across all trials for each subject
    leaveT = readtable(sprintf('%s_trialbytrial.csv', study{s}));
    if strcmp(study{s},'contrerashuerta')
        leaveT = leaveT(leaveT.ben == 1,:); % exclude other condition
    end

    subjSD = zeros(1, numel(unique(leaveT.sub)));
    for iS = unique(leaveT.sub)'
        subjSD(iS) = std(leaveT.leaveT(leaveT.sub == iS,:));
    end

    % get sd leave times across all trials for each simulation
    load(sprintf('modelLT_all_%s_M27', study{s}), '-mat')

    simSD = simulated_leave_times.sd;

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on
    xlabel('Subject SD_{leave} (s)');
    ylabel('Model SD_{leave} (s)');

    scatter(subjSD, simSD, 50, "MarkerFaceColor",color.general, "MarkerEdgeColor",'w')
    if strcmp(study{s},'kane') % set manually for Kane since subjects SD so similar
        xlim([3 4]), ylim([3 4])
    else
        xlim([round(min(subjSD))-1 round(max(subjSD))+1]), ylim([round(min(simSD))-1 round(max(simSD))+1])
    end

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig4_%s_SD_per_participant', study{s})],'-dsvg')
        print([overleaf_path, 'fig4/', sprintf('fig4_%s_SD_per_participant', study{s})],'-dsvg')
    end

    % statistics
    % data not normal for CH and LH
    %     figure, hist(subjSD); [h,p] = swtest(subjSD) 

    %     figure, hist(simSD); [h,p] = swtest(simSD)

    if strcmp(study{s}, 'Kane')
        [r p] = corr(subjSD',simSD', 'type', 'Pearson')
    else 
        [r p] = corr(subjSD',simSD', 'type', 'Spearman')
    end
end



