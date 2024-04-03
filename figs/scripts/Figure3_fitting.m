%%% Figure 3 foraging paper 
% Emma Scholey 7 Aug 2023

clearvars; close all

addpath('./plotFunctions')
addpath('../../model/helperFunctions/')
addpath('../../data/fitting_data_240117/')
addpath('../../model/helperFunctions/')
addpath('../../data/experiment_data/')
addpath('../../data/simulation_data/')
addpath('../../data/analytical_data/')

run figure_properties_foraging.m

study = {'leheron', 'contrerashuerta', 'kane'};

save_figs = 1; % whether to save figures

%% Panel: BIC comparison for model

model = [1, 26:29]; % model numbers to compare
modelNames = {'vary \beta', 'vary \beta vary c', 'vary \beta fix c', 'fix \beta vary c', 'fix \beta fix c'};

nModels = size(model,2);

models_AIC = zeros(nModels, 1);
models_BIC = zeros(nModels, 1);

for s = 1:numel(study)

    for m = 1:nModels
        load(sprintf('fitting_results_M%d_%s', model(m),study{s}), '-mat', 'BIC');
        ppts_BIC(:,:,m) = BIC;
        models_BIC(m,:) = sum(BIC);
    end

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    bar(models_BIC,'FaceColor',color.general, 'EdgeColor', color.general); hold on

    % find best model and highlight
    minBIC = min(models_BIC);
    bestModel = find(models_BIC == minBIC);
    bar(bestModel,minBIC, 'FaceColor',color.highlight, 'EdgeColor', color.highlight)

    ylabel('BIC (sum)')
    set(gca,'XTickLabel',modelNames, 'XTickLabelRotation',50, 'YLim', [min(models_BIC)-200,max(models_BIC)+200])

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig3_%s_BIC', study{s})],'-dsvg')
        print([overleaf_path, 'fig3/', sprintf('fig3_%s_BIC', study{s})],'-dsvg')
    end
    % compute and plot posterior probabilities for each subject
    posteriorProbabilities(:,:,1) = BICposterior(squeeze(ppts_BIC(:,:,[3,1]))); % difference between winning model, and each other model
    posteriorProbabilities(:,:,2) = BICposterior(squeeze(ppts_BIC(:,:,[3,2])));
    posteriorProbabilities(:,:,3) = BICposterior(squeeze(ppts_BIC(:,:,[3,4])));
    posteriorProbabilities(:,:,4) = BICposterior(squeeze(ppts_BIC(:,:,[3,5])));

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.vertical);
    tl = tiledlayout(4, 1);
    for ii = 1:4
        nexttile;
        posteriors = posteriorProbabilities(:,1,ii) - posteriorProbabilities(:,2,ii);
        bar(1:size(posteriorProbabilities,1), sort(posteriors), 'FaceColor',color.general, 'EdgeColor', color.general)
        ylabel('p(model) difference'), ylim([-1.2, 1.2]), yticks([-1:0.5:1])
        xticks([0,numel(posteriors)]), xlabel('Subjects')
    end

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig3_%s_subject_posteriors', study{s})],'-dsvg')
        print([overleaf_path, 'fig3/', sprintf('fig3_%s_subject_posteriors', study{s})],'-dsvg')
    end
    clear ppts_BIC models_BIC posteriorProbabilities
end

%% Panels: simulated vs empirical mean leaving times
% for winning model, M27
for s = 1:numel(study)
    task = buildTask(study{s});

    load(sprintf('fitting_results_M27_%s',study{s}), '-mat');
    load(sprintf('subjectLT_%s', study{s}), '-mat');
    load(sprintf('modelLT_%s_M27', study{s}), '-mat')

    % summarise data
    modelMean = mean(simulated_leave_times.mean,3); modelSEM = std(simulated_leave_times.mean, [], 3)./sqrt(size(simulated_leave_times.mean,3));
    subjectMean = mean(subject_leave_times.mean,3); subjectSEM = std(subject_leave_times.mean, [], 3)./sqrt(size(subject_leave_times.mean,3));

    % Mean panel
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.rectangle)
    tl = tiledlayout(1,2, 'TileSpacing', 'tight');
    nexttile
    % plot the actual experimental data
    errorbar(subjectMean(1,:),subjectSEM(1,:),lines.exp, 'Color',color.rich ,'LineWidth', widths.plot); hold on
    errorbar(subjectMean(2,:),subjectSEM(2,:),lines.exp, 'Color',color.poor,'LineWidth', widths.plot)

    set(gca,'XLim',[0, task.nPatch+1],'XTick',1:task.nPatch,'XTickLabel',task.patchNames, 'YLim', [0 round(max(modelMean, [], 'all')+5,-1)])
    ylabel('Patch leaving time (s)')

    % plot the simulated data on another panel
    nexttile
    errorbar(modelMean(1,:),modelSEM(1,:),lines.model, 'Color',color.rich ,'LineWidth', widths.plot), hold on
    errorbar(modelMean(2,:),modelSEM(2,:),lines.model, 'Color',color.poor,'LineWidth', widths.plot)
    set(gca,'XLim',[0, task.nPatch+1],'XTick',1:task.nPatch,'XTickLabel',task.patchNames, 'YLim', [0 round(max(modelMean, [], 'all')+5,-1)])
    set(gca, 'YTick', [], 'YColor', 'none');

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig3_%s_simulated_mean_LT', study{s})],'-dsvg')
        print([overleaf_path, 'fig3/', sprintf('fig3_%s_simulated_mean_LT', study{s})],'-dsvg')
    end

    % plot on one panel to show overlap, and add MVT
     % Mean panel
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square)
    nexttile
    % plot the actual experimental data
    errorbar(subjectMean(1,:),subjectSEM(1,:),lines.exp, 'Color',color.rich ,'LineWidth', widths.plot); hold on
    errorbar(subjectMean(2,:),subjectSEM(2,:),lines.exp, 'Color',color.poor,'LineWidth', widths.plot)

    set(gca,'XLim',[0, task.nPatch+1],'XTick',1:task.nPatch,'XTickLabel',task.patchNames, 'YLim', [0 round(max(modelMean, [], 'all')+5,-1)])
    ylabel('Patch leaving time (s)')

    errorbar(modelMean(1,:),modelSEM(1,:),lines.model, 'Color',color.rich ,'LineWidth', widths.plot)
    errorbar(modelMean(2,:),modelSEM(2,:),lines.model, 'Color',color.poor,'LineWidth', widths.plot)

    % also plot MVT
    plot(1:task.nPatch,task.optLT(1,:),':','Color',color.rich,'LineWidth',widths.plot)
    plot(1:task.nPatch,task.optLT(2,:),':','Color',color.poor,'LineWidth',widths.plot)

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig3_%s_simulated_mean_LT_overlap', study{s})],'-dsvg')
        print([overleaf_path, 'fig3/', sprintf('fig3_%s_simulated_mean_LT_overlap', study{s})],'-dsvg')
    end


end

%% Panel: per-participant fit (simulated vs actual) 
% for winning model, M27
for s = 1:numel(study)
    task = buildTask(study{s});

    % get mean leave times across all trials for each subject
    leaveT = readtable(sprintf('%s_trialbytrial.csv', study{s}));
    if strcmp(study{s},'contrerashuerta')
        leaveT = leaveT(leaveT.ben == 1,:); % exclude other condition
    end

    subjMean = zeros(1, numel(unique(leaveT.sub)));
    for iS = unique(leaveT.sub)'
        subjMean(iS) = mean(leaveT.leaveT(leaveT.sub == iS,:));
    end

    % get mean leave times across all trials for each simulation
    load(sprintf('modelLT_all_%s_M27', study{s}), '-mat')

    simMean = simulated_leave_times.mean;

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on
    xlabel('Subject leaving times (s)');
    ylabel('Model leaving times (s)');

    scatter(subjMean, simMean, 50, "MarkerFaceColor",color.general, "MarkerEdgeColor",'w')
    if strcmp(study{s},'kane') % set manually for Kane since subjects SD so similar
        xlim([5 10]), ylim([5 10])
    else
        xlim([0 round(max(subjMean))+3]), ylim([0 round(max(subjMean))+3])
    end

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig3_%s_per_participant', study{s})],'-dsvg')
        print([overleaf_path, 'fig3/', sprintf('fig3_%s_per_participant', study{s})],'-dsvg')
    end

    % statistics
    % check distributions
%     figure, hist(subjMean)
%     figure, hist(simMean)

    [r p] = corr(subjMean',simMean')
end


    %% Supplementary Panel: beta fit x subject LT 
% for winning model, M27
for s = 1:numel(study)
    task = buildTask(study{s});

    % get mean leave times per environment for each subject
    leaveT = readtable(sprintf('%s_trialbytrial.csv', study{s}));
    if strcmp(study{s},'contrerashuerta')
        leaveT = leaveT(leaveT.ben == 1,:); % exclude other condition
    end

    subjMean = zeros(2, numel(unique(leaveT.sub)));
    for iS = unique(leaveT.sub)'
        for iE = 1:task.nEnviron
            subjMean(iE,iS) = mean(leaveT.leaveT(leaveT.sub == iS & leaveT.env == iE,:));
        end
    end

    % get beta fits for each subject
    load(sprintf('fitting_results_M27_%s',study{s}), '-mat', 'minNLLFitParams');

    % get subject LT per environment

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on
    xlabel('Subject leaving times (s)');
    ylabel('Estimated beta fit');

    scatter(subjMean(1,:), minNLLFitParams.beta_rich, 50, "MarkerFaceColor",color.rich, "MarkerEdgeColor",'w','MarkerFaceAlpha', 0.5,'MarkerEdgeAlpha',0.5), hold on
    scatter(subjMean(2,:), minNLLFitParams.beta_poor, 50, "MarkerFaceColor",color.poor, "MarkerEdgeColor",'w','MarkerFaceAlpha', 0.5,'MarkerEdgeAlpha',0.5)

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig3_%s_beta_LT_supplement', study{s})],'-dsvg')
        print([overleaf_path, 'fig3/', sprintf('fig3_%s_beta_LT_supplement', study{s})],'-dsvg')
    end

    % statistics
        [h,p] = swtest([minNLLFitParams.beta_rich; minNLLFitParams.beta_poor]) % data not normal
% figure, hist(log([minNLLFitParams.beta_rich; minNLLFitParams.beta_poor]))
        [r p] = corr([subjMean(1,:),subjMean(2,:)]',[minNLLFitParams.beta_rich; minNLLFitParams.beta_poor], 'type', 'Spearman')

end


%% Panel: distribution of beta fit with LT for each environment
% for winning model, M27
for s = 1:numel(study)

    load(sprintf('fitting_results_M27_%s',study{s}), '-mat', 'minNLLFitParams');

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on

    scatter(1, minNLLFitParams.beta_rich, 70, 'MarkerFaceColor',color.rich, 'MarkerEdgeColor','w'); hold on
    scatter(2, minNLLFitParams.beta_poor, 70, 'MarkerFaceColor',color.poor, 'MarkerEdgeColor','w')

    % connect points
    plot([1,2], [minNLLFitParams.beta_rich,minNLLFitParams.beta_poor], 'Color', [0.7 0.7 0.7])

    errorbar(1, median(minNLLFitParams.beta_rich), std(minNLLFitParams.beta_rich)./sqrt(numel(minNLLFitParams.beta_rich)), "_", "CapSize", 0, "LineWidth", widths.plot, "Color", 'k')
    errorbar(2, median(minNLLFitParams.beta_poor), std(minNLLFitParams.beta_poor)./sqrt(numel(minNLLFitParams.beta_rich)), "_", "CapSize", 0, "LineWidth", widths.plot, "Color", 'k')

    ylabel('\beta (higher = exploit)')
    if strcmp(study{s},'kane')
      set(gca,'XLim', [0 3],'XTick', [1:2], 'XTickLabel', {'Rich', 'Poor'})
    else
        set(gca,'XLim', [0 3],'XTick', [1:2], 'XTickLabel', {'Rich', 'Poor'})
    end

    % statistics
%     [h,p] = swtest(minNLLFitParams.beta_rich) % data not normal even if we log transform
%     figure, hist(minNLLFitParams.beta_rich)
%     [h,p] = swtest(minNLLFitParams.beta_poor) % data normal if we log transform
%         figure, hist(minNLLFitParams.beta_poor)
    % use robust test
    [p_val,H0_reject, stats] = signrank(minNLLFitParams.beta_rich, minNLLFitParams.beta_poor, 'method', 'approximate')

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig3_%s_beta_fits', study{s})],'-dsvg')
        print([overleaf_path, 'fig3/', sprintf('fig3_%s_beta_fits', study{s})],'-dsvg')
    end
end


%% STATISTICS %%%%%%%%%%%%%%%%%%%%
% % statistics - mean LT - confirms le heron results
% rich_mean = squeeze(subject_leave_times.mean(1,:,:))';
% poor_mean = squeeze(subject_leave_times.mean(2,:,:))';
% 
% mean_table = [rich_mean;poor_mean];
% 
% %%% test assumptions
% 
% % normality - fine
% figure, hist(rich_mean(:,1))
% % homogeneity of  variance - fine
% vartestn([rich_mean,poor_mean]);
% 
% [p_val, tbl] = anova2(mean_table,size(rich_mean,1));
% 
% % statistics - SD LT - confirms novel SD predictions
% rich_sd = squeeze(subject_leave_times.sd(1,:,:))';
% poor_sd = squeeze(subject_leave_times.sd(2,:,:))';
% 
% sd_table = [rich_sd;poor_sd];
% 
% %normality - not fine
% figure, hist(rich_sd(:,3))
% figure, hist(poor_sd(:,3))
% % homogeneity of  variance - fine
% vartestn([rich_sd,poor_sd]);
% 
% [p_val, tbl] = anova2(sd_table,size(rich_sd,1));
% 
% % statistics - SD LT SIMULATIONS
% rich_sd = squeeze(subject_leave_times_simulated.sd(1,:,:))';
% poor_sd = squeeze(subject_leave_times_simulated.sd(2,:,:))';
% 
% sd_table = [rich_sd;poor_sd];
% 
% [p_val, tbl] = anova2(sd_table,size(rich_sd,1));
