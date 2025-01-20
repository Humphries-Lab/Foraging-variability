%%% Figure 3 foraging paper
% Emma Scholey 7 Aug 2023

clearvars; close all

addpath('./plotFunctions')
addpath('../../model/helperFunctions/')
addpath('../../data/fitting_data/')
addpath('../../model/helperFunctions/')
addpath(genpath('../../data/experiment_data/'))
addpath('../../data/simulation_data/')
addpath('../../data/analytical_data/')

run figure_properties_foraging.m

study = {'leheron', 'contrerashuerta', 'kane'};

save_figs = 1; % whether to save figures

%% Panel: BIC comparison for model

model = [1:5]; % model numbers to compare
modelNames = {'vary \beta', 'vary \beta vary c', 'vary \beta fix c', 'fix \beta vary c', 'fix \beta fix c'};

nModels = size(model,2);

models_AIC = zeros(nModels, 1);
models_BIC = zeros(nModels, 1);

for s = 1:numel(study)

    for m = 1:nModels
        load(sprintf('fitting_results_M%d_%s', model(m),study{s}), '-mat', 'BIC');
        ppts_BIC(:,m) = BIC;
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
        print([export_path, sprintf('fig4_%s_BIC', study{s})],'-dsvg')
    end

    % compute and plot posterior probabilities for each subject
    % posteriorProbabilities = BICposterior(ppts_BIC); % difference between winning model, and each other model
    %
    % figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.rectangle);
    % tl = tiledlayout(1,5);
    % for ii = 1:nModels
    %     nexttile;
    %     histogram(posteriorProbabilities(:,ii), 'FaceColor',color.general, 'EdgeColor', color.general, 'NumBins',5)
    %     xlabel(['p(', modelNames{ii},')']),xticks([0:0.5:1]), xlim([-0.1,1.1])
    %     ylabel('Number of participants'),ylim([0,size(ppts_BIC,1)])
    % end

    % show best model per participant
    subjectBestModelBIC = ppts_BIC == min(ppts_BIC, [], 2);
    sumBestModel = sum(subjectBestModelBIC);

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);

    bar(1:nModels,sumBestModel,'FaceColor',color.general, 'EdgeColor', color.general), hold on
    ylabel('Number of subjects')
    set(gca,'XTickLabel',modelNames, 'XTickLabelRotation',50)
    bar(bestModel,sumBestModel(bestModel), 'FaceColor',color.highlight, 'EdgeColor', color.highlight)

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig4_%s_BIC_subjects', study{s})],'-dsvg')
    end
    clear ppts_BIC models_BIC posteriorProbabilities
end

%% Panels: simulated vs empirical mean leaving times
% for winning model, M3
for s = 1:numel(study)
    task = buildTask(study{s});

    load(sprintf('fitting_results_M3_%s',study{s}), '-mat');
    load(sprintf('subject_LT_%s', study{s}), '-mat');
    load(sprintf('sim_LT_%s_M3', study{s}), '-mat')

    % summarise data
    modelMean = mean(simulated_leave_times.mean,3); modelSEM = std(simulated_leave_times.mean, [], 3)./sqrt(size(simulated_leave_times.mean,3));
    subjectMean = mean(subject_leave_times.mean,3); subjectSEM = std(subject_leave_times.mean, [], 3)./sqrt(size(subject_leave_times.mean,3));

    % plot on one panel to show overlap, and add MVT
    % Mean panel
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square)
    nexttile
    % plot the actual experimental data
    errorbar(subjectMean(1,:),subjectSEM(1,:),lines.exp, 'Color',color.rich ,'LineWidth', widths.plot); hold on
    errorbar(subjectMean(2,:),subjectSEM(2,:),lines.exp, 'Color',color.poor,'LineWidth', widths.plot)

    set(gca,'XLim',[0, task.nPatch+1],'XTick',1:task.nPatch,'XTickLabel',task.patchNames, 'YLim', [0 round(max(modelMean, [], 'all')+5,-1)])
    
    if strcmp(study{s}, 'kane')
        ylabel('Harvests per patch')
    else
        ylabel('Patch leaving time (s)')
    end

    errorbar(modelMean(1,:),modelSEM(1,:),lines.model, 'Color',color.rich ,'LineWidth', widths.plot)
    errorbar(modelMean(2,:),modelSEM(2,:),lines.model, 'Color',color.poor,'LineWidth', widths.plot)

    % also plot MVT
    plot(1:task.nPatch,task.optLT(1,:),lines.mvt,'Color',color.rich,'LineWidth',widths.plot)
    plot(1:task.nPatch,task.optLT(2,:),lines.mvt,'Color',color.poor,'LineWidth',widths.plot)

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig4_%s_simulated_mean_LT_overlap', study{s})],'-dsvg')
    end


    %% Panel: per-participant fit (simulated vs actual)
    % for winning model, M3

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on
    tl = tiledlayout(1,2,'TileSpacing','tight');

    plot_colour = [color.rich;color.poor];

    corr_coord = [0.25,0.65];

    for iE = 1:task.nEnviron
        nexttile, hold on

        tmp_subject = squeeze(subject_leave_times.mean(iE,:,:));
        tmp_model = squeeze(simulated_leave_times.mean(iE,:,:));
        
        if strcmp(study{s},'kane')
            xlim([0 14]), ylim([0 14]),xticks([0:7:14]),yticks([0:7:14])
        elseif strcmp(study{s}, 'contrerashuerta')
            xlim([0 30]), ylim([0 30]), xticks([0:15:30]),yticks([0:15:30])
        elseif strcmp(study{s}, 'leheron')
            xlim([0 60]), ylim([0 60]), xticks([0:30:60]),yticks([0:30:60])
        end

        if iE == 2
            yticklabels('');
        end

        plot(xlim,ylim, ('-k')) % add identity line
        scatter(tmp_subject(:), tmp_model(:), 50, "MarkerFaceColor",plot_colour(iE,:), "MarkerEdgeColor",'w','MarkerFaceAlpha', 0.5,'MarkerEdgeAlpha',0.5)
        
        if iE == 2
            yticklabels('');
        end

        [r p] = corr(tmp_subject(:),tmp_model(:), 'type', 'Pearson')

        str=sprintf('r = %1.2f',r);
        T = textbypos(corr_coord(iE),0.25,str);

    end

    if strcmp(study{s}, 'kane')
        xlabel(tl,'Subject harvests per patch')
        ylabel(tl, 'Model harvests per patch');
    else
        xlabel(tl, 'Subject leaving times (s)');
        ylabel(tl, 'Model leaving times (s)');
    end

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig4_%s_per_participant', study{s})],'-dsvg')
    end
end


%% Supplementary Panel: beta fit x subject LT
% for winning model, M3
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
    load(sprintf('fitting_results_M3_%s',study{s}), '-mat', 'minNLLFitParams');

    % get subject LT per environment

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on
    
    if strcmp(study{s}, 'kane')
        xlabel('Subject harvests per patch')
    else
        xlabel('Subject leaving times (s)');
    end

    ylabel('Estimated beta fit');

    scatter(subjMean(1,:), minNLLFitParams.beta_rich, 50, "MarkerFaceColor",color.rich, "MarkerEdgeColor",'w','MarkerFaceAlpha', 0.5,'MarkerEdgeAlpha',0.5), hold on
    scatter(subjMean(2,:), minNLLFitParams.beta_poor, 50, "MarkerFaceColor",color.poor, "MarkerEdgeColor",'w','MarkerFaceAlpha', 0.5,'MarkerEdgeAlpha',0.5)

    [r p] = corr([subjMean(1,:),subjMean(2,:)]',[minNLLFitParams.beta_rich; minNLLFitParams.beta_poor], 'type', 'Spearman');

    str=sprintf('r_s = %1.2f',r);
    T = textbypos(0.25,0.87,str);

        if p < .001
            sig_val = textbypos(0.25,0.77,'p < .001');
        elseif p >= .001 && p < .01
            sig_val = textbypos(0.25,0.77,sprintf('p = %1.3f', p));
        else
            sig_val = textbypos(0.25,0.77,sprintf('p = %1.2f', p));
        end


    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig4_%s_beta_LT_supplement', study{s})],'-dsvg')
    end

end


%% Panel: distribution of beta fit for each environment
% for winning model, M3
for s = 1:numel(study)

    load(sprintf('fitting_results_M3_%s',study{s}), '-mat', 'minNLLFitParams');

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on

    % connect points
    plot([1,2], [minNLLFitParams.beta_rich,minNLLFitParams.beta_poor], 'Color', [0.8 0.8 0.8]); hold on

    scatter(1, minNLLFitParams.beta_rich, 70, 'MarkerFaceColor',color.rich, 'MarkerEdgeColor','w', 'LineWidth',1)
    scatter(2, minNLLFitParams.beta_poor, 70, 'MarkerFaceColor',color.poor, 'MarkerEdgeColor','w', 'LineWidth',1)

    errorbar(1, mean(minNLLFitParams.beta_rich), std(minNLLFitParams.beta_rich)./sqrt(numel(minNLLFitParams.beta_rich)), "_", "CapSize", 0, "LineWidth", widths.plot, "Color", 'k')
    errorbar(2, mean(minNLLFitParams.beta_poor), std(minNLLFitParams.beta_poor)./sqrt(numel(minNLLFitParams.beta_rich)), "_", "CapSize", 0, "LineWidth", widths.plot, "Color", 'k')

    ylabel('\beta (higher = exploit)')
    set(gca,'XLim', [0 3],'XTick', [1:2], 'XTickLabel', {'Rich', 'Poor'})

    % statistics
    %figure, qqplot(minNLLFitParams.beta_poor - minNLLFitParams.beta_rich)
     [~,p_sw] = swtest(minNLLFitParams.beta_poor - minNLLFitParams.beta_rich);

    if p_sw < .05 | strcmp(study{s}, 'kane')
        [p, ~, stats] = signrank(minNLLFitParams.beta_poor,minNLLFitParams.beta_rich, 'method', 'exact');
    else
        [~, p, ~, stats] = ttest(minNLLFitParams.beta_poor,minNLLFitParams.beta_rich)
    end

        if p < .001
            sig_val = textbypos(0.5, 0.8,'p < .001');
        elseif p >= .001 && p < .01
            sig_val = textbypos(0.5, 0.8,sprintf('p = %1.3f', p));
        else
            sig_val = textbypos(0.5, 0.8,sprintf('p = %1.2f', p));
        end


    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig4_%s_beta_fits', study{s})],'-dsvg')
    end

end


