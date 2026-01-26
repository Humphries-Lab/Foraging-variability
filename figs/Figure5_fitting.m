%%% Figure 5 - model fitting and comparison foraging paper
% Emma Scholey 7 Aug 2023
% last updated 6 October 2025

clearvars; close all

addpath('./plotFunctions')
addpath('../model/helperFunctions/')
addpath(genpath('../data/'))

run figure_properties_foraging.m

study = {'leheron', 'contrerashuerta', 'kane'};

save_figs = 1; % whether to save figures

%% Panel: BIC comparison for model

model = [1:5]; % model numbers to compare

modelNames = {'vary B', 'vary B vary c', 'vary B fix c', 'fix B vary c', 'fix B fix c'};

nModels = size(model,2);

models_AIC = zeros(nModels, 1);
models_BIC = zeros(nModels, 1);

for s = 1:numel(study)

    for m = 1:nModels
        load(sprintf('fitting_results_M%d_%s', model(m),study{s}), '-mat', 'BIC');
        ppts_BIC(:,m) = BIC;
        models_BIC(m,:) = sum(BIC);
    end

     % find best model and highlight
    minBIC = min(models_BIC);
    bestModel = find(models_BIC == minBIC);

    diff_BIC = models_BIC - minBIC;
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    bar(diff_BIC,'FaceColor',color.general, 'EdgeColor', color.general); hold on

   
    %bar(bestModel,minBIC, 'FaceColor',color.highlight, 'EdgeColor', color.highlight)

    ylabel('Diff. from lowest BIC')
    set(gca,'XTickLabel',modelNames, 'XTickLabelRotation',50)

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig5_%s_BIC', study{s})],'-dsvg')
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
        print([export_path, sprintf('fig5_%s_BIC_subjects', study{s})],'-dsvg')
    end
    clear ppts_BIC models_BIC posteriorProbabilities
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

    % figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on
    % 
    % if strcmp(study{s}, 'kane')
    %     xlabel('Subject harvests per patch')
    % else
    %     xlabel('Subject leaving times (s)');
    % end
    % 
    % ylabel('Rew. weight, B');
    % 
    % scatter(subjMean(1,:), minNLLFitParams.beta_rich, 50, "MarkerFaceColor",color.rich, "MarkerEdgeColor",'w','MarkerFaceAlpha', 0.5,'MarkerEdgeAlpha',0.5), hold on
    % scatter(subjMean(2,:), minNLLFitParams.beta_poor, 50, "MarkerFaceColor",color.poor, "MarkerEdgeColor",'w','MarkerFaceAlpha', 0.5,'MarkerEdgeAlpha',0.5)
    % 
    % [r p] = corr([subjMean(1,:),subjMean(2,:)]',[minNLLFitParams.beta_rich; minNLLFitParams.beta_poor], 'type', 'Spearman');
    % 
    % str=sprintf('r_s = %1.2f',r);
    % T = textbypos(0.25,0.87,str);
    % 
    % if p < .001
    %     sig_val = textbypos(0.25,0.77,'p < .001');
    % elseif p >= .001 && p < .01
    %     sig_val = textbypos(0.25,0.77,sprintf('p = %1.3f', p));
    % else
    %     sig_val = textbypos(0.25,0.77,sprintf('p = %1.2f', p));
    % end
    % 
    % 
    % FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    % if save_figs == 1
    %     print([export_path, sprintf('fig5_%s_beta_LT', study{s})],'-dsvg')
    % end

    %% Panel: beta difference versus environment difference correlation
    % for winning model, M3

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on

    beta_diff = minNLLFitParams.beta_poor - minNLLFitParams.beta_rich;
    env_diff = subjMean(2,:) - subjMean(1,:);

    scatter(beta_diff, env_diff, 50, 'MarkerFaceColor',color.general, "MarkerEdgeColor",'w','MarkerFaceAlpha', 0.5,'MarkerEdgeAlpha',0.5)

    xlabel('\Delta B: poor - rich')
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

    if strcmp(study{s}, 'kane')
        ylim([0 inf])
    end

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig5_%s_env_beta_diffs', study{s})],'-dsvg')
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

    ylabel('Bias, c');

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
        print([export_path, sprintf('fig5_%s_c_LT', study{s})],'-dsvg')
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

    ylabel('Rew. weight, B')
    set(gca,'XLim', [0 3],'XTick', [1:2], 'XTickLabel', {'Rich', 'Poor'})

    % statistics
    %figure, qqplot(minNLLFitParams.beta_poor - minNLLFitParams.beta_rich)
    [~,p_sw] = swtest(minNLLFitParams.beta_poor - minNLLFitParams.beta_rich);

    if p_sw < .05 | strcmp(study{s}, 'kane')
        [p, ~, stats] = signrank(minNLLFitParams.beta_poor,minNLLFitParams.beta_rich, 'method', 'exact');
        meanEffectSize(minNLLFitParams.beta_poor,minNLLFitParams.beta_rich,Paired=true, Effect='robustcohen')

    else
        [~, p, ~, stats] = ttest(minNLLFitParams.beta_poor,minNLLFitParams.beta_rich)
        meanEffectSize(minNLLFitParams.beta_poor,minNLLFitParams.beta_rich,Paired=true, Effect='cohen')

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
        print([export_path, sprintf('fig5_%s_beta_fits', study{s})],'-dsvg')
    end

end


%% Panel: model identifiability

for s = 1:numel(study)

    load(sprintf('XP_model_identifiability_%s', study{s}), 'exceedance_probabilities', 'protected_XP')
    modelNames = {'vary B', 'vary B vary c', 'vary B fix c', 'fix B vary c', 'fix B fix c'};

    exceedance_probabilities = round(exceedance_probabilities,2);

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square)
    h = heatmap(exceedance_probabilities,'MissingDataColor','w', 'ColorMap', sky, 'GridVisible', 'off','CellLabelColor','none');
    labels = modelNames;
    h.XDisplayLabels = labels;
    h.YDisplayLabels = labels;
    h.XLabel = 'estimated';
    h.YLabel = 'simulated';
    q = struct(h);
    q.XAxis.TickLabelRotation = 60;
    set(gca,'FontSize',fontsize)
    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig5_%s_model_identifiability', study{s})],'-dsvg')
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
    h.XDisplayLabels = {'B_{rich}', 'B_{poor}', 'c'};
    h.YDisplayLabels = {'B_{rich}', 'B_{poor}', 'c'};
    h.XLabel = 'estimated';
    h.YLabel = 'simulated';
    %set(gca,'FontSize',fontsize)

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig5_%s_param_recovery', study{s})],'-dsvg')
    end

end


