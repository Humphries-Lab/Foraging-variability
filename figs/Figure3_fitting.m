%%% Figure 3 foraging paper
% Emma Scholey 7 Aug 2023
% last updated 11 June 2025

clearvars; close all

addpath('./plotFunctions')
addpath('../model/helperFunctions/')
addpath(genpath('../data/'))

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
        print([export_path, sprintf('fig3_%s_BIC', study{s})],'-dsvg')
    end

    % show best model per participant
    subjectBestModelBIC = ppts_BIC == min(ppts_BIC, [], 2);
    sumBestModel = sum(subjectBestModelBIC);
    bestModel = find(sumBestModel == max(sumBestModel));

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);

    bar(1:nModels,sumBestModel,'FaceColor',color.general, 'EdgeColor', color.general), hold on
    ylabel('Number of subjects')
    set(gca,'XTickLabel',modelNames, 'XTickLabelRotation',50)
    bar(bestModel,sumBestModel(bestModel), 'FaceColor',color.highlight, 'EdgeColor', color.highlight)

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig3_%s_BIC_subjects', study{s})],'-dsvg')
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
        print([export_path, sprintf('fig3_%s_simulated_mean_LT_overlap', study{s})],'-dsvg')
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
        print([export_path, sprintf('fig3_%s_per_participant', study{s})],'-dsvg')
    end
end


%% Panel: beta fit x subject LT
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
        print([export_path, sprintf('fig3_%s_beta_LT_supplement', study{s})],'-dsvg')
    end

    %% Panel: beta difference versus environment difference correlation
    % for winning model, M3

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on

    beta_diff = minNLLFitParams.beta_poor - minNLLFitParams.beta_rich;
    env_diff = subjMean(2,:) - subjMean(1,:);

    scatter(beta_diff, env_diff, 50, 'MarkerFaceColor',color.general, "MarkerEdgeColor",'w','MarkerFaceAlpha', 0.5,'MarkerEdgeAlpha',0.5)

    xlabel('\Delta \beta: poor - rich')
    if strcmp(study{s}, 'kane')
        ylabel('\Delta harvests: poor - rich')
    else
        ylabel('\Delta leaving times: poor - rich')
    end

    [r p] = corr(beta_diff, env_diff', 'type', 'Spearman');

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
        print([export_path, sprintf('fig3_%s_env_beta_diffs', study{s})],'-dsvg')
    end

end

%% Panel: C fit x subject LT
% for winning model, M3
for s = 1:numel(study)
    task = buildTask(study{s});

    % get mean leave times per environment for each subject
    leaveT = readtable(sprintf('%s_trialbytrial.csv', study{s}));
    if strcmp(study{s},'contrerashuerta')
        leaveT = leaveT(leaveT.ben == 1,:); % exclude other condition
    end

    subjMean = zeros(1, numel(unique(leaveT.sub)));
    for iS = unique(leaveT.sub)'
        subjMean(iS) = mean(leaveT.leaveT(leaveT.sub == iS,:));
    end

    % get beta fits for each subject
    load(sprintf('fitting_results_M3_%s',study{s}), '-mat', 'minNLLFitParams');

    % get subject LT per environment

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on
    l = yline(0, '--', 'LineWidth',1, 'Color','k');

    if strcmp(study{s}, 'kane')
        xlabel('Subject harvests per patch')
    else
        xlabel('Subject leaving times (s)');
    end

    ylabel('Estimated c fit');

    scatter(subjMean, minNLLFitParams.bias, 50, "MarkerFaceColor",color.general, "MarkerEdgeColor",'w','MarkerFaceAlpha', 0.5,'MarkerEdgeAlpha',0.5)

    % [r p] = corr(subjMean',minNLLFitParams.bias, 'type', 'Spearman');
    % 
    % str=sprintf('r_s = %1.2f',r);
    % T = textbypos(0.6,0.4,str);
    % 
    % if p < .001
    %     sig_val = textbypos(0.6,0.32,'p < .001');
    % elseif p >= .001 && p < .01
    %     sig_val = textbypos(0.6,0.32,sprintf('p = %1.3f', p));
    % else
    %     sig_val = textbypos(0.6,0.32,sprintf('p = %1.2f', p));
    % end


    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig3_%s_c_LT_supplement', study{s})],'-dsvg')
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

    ylabel('\beta')
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
        print([export_path, sprintf('fig3_%s_beta_fits', study{s})],'-dsvg')
    end

end


%% Panel: model identifiability

for s = 1:numel(study)

    load(sprintf('XP_model_identifiability_%s', study{s}), 'exceedance_probabilities', 'protected_XP', 'modelNames')

    exceedance_probabilities = round(exceedance_probabilities,2);

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square)
    h = heatmap(exceedance_probabilities,'MissingDataColor','w', 'ColorMap', sky, 'GridVisible', 'off','CellLabelColor','none');
    labels = modelNames;
    h.XDisplayLabels = labels;
    h.YDisplayLabels = labels;
    h.XLabel = 'estimated';
    h.YLabel = 'simulated';
    q = struct(h);
    q.XAxis.TickLabelRotation = 50;
    set(gca,'FontSize',fontsize)
    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig3_%s_model_identifiability', study{s})],'-dsvg')
    end

end


%% Panel: parameter recovery
% for winning model only

modelNum = [3];

for s = 1:numel(study)

    load(sprintf('M%d_parameter_recovery_%s', modelNum, study{s}), 'params')

    r = corr(params.simulated, params.estimated, type="Spearman");
    r = round(r,2);
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square)
    h = heatmap(r,'MissingDataColor','w', 'ColorMap', sky, 'GridVisible', 'off');
    labels = params.names;
    h.XDisplayLabels = labels;
    h.YDisplayLabels = labels;
    h.XLabel = 'estimated';
    h.YLabel = 'simulated';
    %set(gca,'FontSize',fontsize)

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig3_%s_param_recovery', study{s})],'-dsvg')
    end

end


%% Panel: expected leaving times

panelData = {
    'rangebeta';...
    'rangebias';...
    };

for s = 1:numel(study)

    task = buildTask(study{s});

    for p = 1:numel(panelData)

        % Panel: expected leave times for range of exploit parameter
        load(sprintf('expectedLT_%s_%s.mat', study{s}, panelData{p}));
        load(sprintf('subject_LT_%s.mat', study{s}));

        E_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
        colororder(color.patch);

        if strcmp(panelData{p},'rangebias')
            signedlog = @(x) sign(x) .* log10(1 + abs(x));
            data.explore = signedlog(data.explore);
            plot(data.explore, data.E_leave, 'LineWidth', widths.plot); hold on

            xlabel('{\it c}:')
            xticks = [-60 -10 -1 0 1 3];  % Choose meaningful bias values
            set(gca, 'XTick', signedlog(xticks))
            set(gca, 'XTickLabel', arrayfun(@num2str, xticks, 'UniformOutput', false))

        elseif strcmp(panelData{p},'rangebeta')
            semilogx(data.explore,data.E_leave,'LineWidth',widths.plot); hold on

            xlabel('\beta:')
            xlim([10^-1.5 10^0])
        end

        if strcmp(study{s}, 'kane')
            ylabel('Expected harvests per patch')
        else
            ylabel('Expected leaving time (s)')
        end


        FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
        % beta prediction for MVT leaving times, based on mid patch

        ixRich_MVT = find(data.E_leave(:,2) >= task.optLT(1,2),1,"first");
        ixPoor_MVT = find(data.E_leave(:,2) >= task.optLT(2,2),1,"first");
        parameter_rich = data.explore(ixRich_MVT);
        parameter_poor = data.explore(ixPoor_MVT);

        line([parameter_rich parameter_rich],[0 E_fig.Children.YLim(2)],'Color',[color.rich,0.8], 'LineWidth',widths.plot, 'LineStyle', lines.mvt)
        line([parameter_poor parameter_poor],[0 E_fig.Children.YLim(2)],'Color',[color.poor,0.8], 'LineWidth',widths.plot, 'LineStyle', lines.mvt)

        mean_subj_LT = mean(subject_leave_times.mean,3);
        % beta prediction for subject leaving times, based on mid patch
        ixRich_subj = find(data.E_leave(:,2) >= mean_subj_LT(1,2),1,"first");
        ixPoor_subj = find(data.E_leave(:,2) >= mean_subj_LT(2,2),1,"first");
        parameter_rich = data.explore(ixRich_subj);
        parameter_poor = data.explore(ixPoor_subj);

        line([parameter_rich parameter_rich],[0 E_fig.Children.YLim(2)],'Color',[color.rich,0.8], 'LineWidth',widths.plot)
        line([parameter_poor parameter_poor],[0 E_fig.Children.YLim(2)],'Color',[color.poor,0.8], 'LineWidth',widths.plot)

        print([export_path 'fig3_E_leave_' study{s} '_' panelData{p}],'-dsvg')

        % Panel: patch leaving time predictions for model

        figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.small_panel);

        plot(1:task.nPatch,data.E_leave(ixRich_MVT,:),lines.mvt,'Color',color.rich,'LineWidth',widths.plot); hold on
        plot(1:task.nPatch,data.E_leave(ixPoor_MVT,:),lines.mvt,'Color',color.poor,'LineWidth',widths.plot)
        plot(1:task.nPatch,data.E_leave(ixRich_subj,:),'Color',color.rich,'LineWidth',widths.plot); hold on
        plot(1:task.nPatch,data.E_leave(ixPoor_subj,:),'Color',color.poor,'LineWidth',widths.plot)

        xlabel('Patch yield')
        if strcmp(study{s}, 'kane')
            ylabel('Harvests per patch')
        else
            ylabel('Leaving time (s)')
        end

        set(gca,'XTick',1:task.nPatch,'XTickLabel',task.patchNames,'XLim',[0 task.nPatch+1],'YLim',[0 round(max(max(mean_subj_LT)))+5])
        if strcmp(study{s}, 'kane')
            ylim([0,14])
        end
        FormatFig_For_Export(gcf,fontsize-1,fontname,widths.axis)
        print([export_path 'fig3_patch_times_' study{s} '_' panelData{p}],'-dsvg')

        clear task.optLT
    end
end

