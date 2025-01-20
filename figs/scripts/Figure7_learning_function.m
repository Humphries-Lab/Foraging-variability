%%% Figure 6 foraging paper
% Emma Scholey 9 April 2024

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

model = [8 11 14 17]; % compare MVT and Block models
modelNames = {'learn \beta M1', 'learn \beta M2', 'learn \beta M3','learn \beta M4'};

nModels = size(model,2);

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
    bestModelIndex = find(models_BIC == minBIC);
    bestModel = model(bestModelIndex);

    bar(bestModelIndex,minBIC, 'FaceColor',color.highlight, 'EdgeColor', color.highlight)

    ylabel('BIC (sum)')
    set(gca, 'XLim', [0,nModels+1], 'XTickLabel',modelNames, 'XTickLabelRotation',50, 'YLim', [min(models_BIC)-200,max(models_BIC)+200])

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig7_%s_BIC', study{s})],'-dsvg')
    end

    % show best model per participant
    subjectBestModelBIC = ppts_BIC == min(ppts_BIC, [], 2);
    sumBestModel = sum(subjectBestModelBIC);

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);

    bar(1:nModels,sumBestModel,'FaceColor',color.general, 'EdgeColor', color.general), hold on
    ylabel('Number of subjects')
    set(gca,'XTickLabel',modelNames, 'XTickLabelRotation',50)

    sumBestModelIndex = find(sumBestModel == max(sumBestModel));
    bar(sumBestModelIndex, sumBestModel(sumBestModelIndex),'FaceColor',color.highlight, 'EdgeColor', color.highlight)

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig7_%s_BIC_subjects', study{s})],'-dsvg')
    end

    clear ppts_BIC

    %% Panels: simulated vs empirical mean leaving times
    % for winning model
    bestModel = 17; % to keep simple, just do next plots for the 2 parameter model
    task = buildTask(study{s});

    load(sprintf('fitting_results_M%d_%s',bestModel,study{s}), '-mat');
    load(sprintf('subject_LT_%s', study{s}), '-mat');
    load(sprintf('sim_LT_%s_M%d', study{s}, bestModel), '-mat')

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

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig7_%s_simulated_mean_LT_overlap', study{s})],'-dsvg')
    end

    %% Panels: simulated vs empirical SD leaving times
    % for winning model

    % summarise data
    modelSD = mean(simulated_leave_times.sd,3); modelSD_SEM = std(simulated_leave_times.sd, [], 3)./sqrt(size(simulated_leave_times.sd,3));
    subjectSD = mean(subject_leave_times.sd,3); subjectSD_SEM = std(subject_leave_times.sd, [], 3)./sqrt(size(subject_leave_times.sd,3));

    % plot on one panel to show overlap, and add MVT
    % Mean panel
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square)
    nexttile
    % plot the actual experimental data
    errorbar(subjectSD(1,:),subjectSD_SEM(1,:),lines.exp, 'Color',color.rich ,'LineWidth', widths.plot); hold on
    errorbar(subjectSD(2,:),subjectSD_SEM(2,:),lines.exp, 'Color',color.poor,'LineWidth', widths.plot)

    set(gca,'XLim',[0, task.nPatch+1],'XTick',1:task.nPatch,'XTickLabel',task.patchNames)
   if strcmp(study{s}, 'kane')
        ylabel('SD_{leave} (harvests)');
    else
        ylabel('SD_{leave} (s)');
   end

    if strcmp(study{s},'kane') % set manually for Kane since subjects SD so similar
        ylim([1 3])
    elseif strcmp(study{s}, 'contrerashuerta')
        ylim([1 4])
    elseif strcmp(study{s}, 'leheron')
        ylim([2 6])
    end

    errorbar(modelSD(1,:),modelSD_SEM(1,:),lines.model, 'Color',color.rich ,'LineWidth', widths.plot)
    errorbar(modelSD(2,:),modelSD_SEM(2,:),lines.model, 'Color',color.poor,'LineWidth', widths.plot)

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig7_%s_simulated_SD_LT_overlap', study{s})],'-dsvg')
    end

    %% Panel: per-participant fit (simulated vs actual) for MEANS
    % for winning model

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on
    tl = tiledlayout(1,2,'TileSpacing','tight');

    plot_colour = [color.rich;color.poor];

    corr_coord = [0.3,0.67];
    if strcmp(study{s}, 'kane')
        xlabel(tl,'Subject harvests per patch')
        ylabel(tl, 'Model harvests per patch');
    else
        xlabel(tl, 'Subject leaving times (s)');
        ylabel(tl, 'Model leaving times (s)');
    end
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
        T = textbypos(corr_coord(iE),0.35,str);

    end



    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig7_%s_per_participant', study{s})],'-dsvg')
    end

    %% Panel: per-participant fit (simulated vs actual) for VARIABILITY

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on
    tl = tiledlayout(1,2,'TileSpacing','tight');

    plot_colour = [color.rich;color.poor];

    corr_coord = [0.15,0.62];

    for iE = 1:task.nEnviron
        nexttile, hold on

        tmp_subject = squeeze(subject_leave_times.sd(iE,:,:));
        tmp_model = squeeze(simulated_leave_times.sd(iE,:,:));

        if strcmp(study{s},'kane')
            xlim([1 3]), ylim([1 3]), xticks([1:1:3]),yticks([1:1:3])
        elseif strcmp(study{s}, 'contrerashuerta')
            xlim([0 8]), ylim([0 8]), xticks([0:4:8]),yticks([0:4:8])
        elseif strcmp(study{s}, 'leheron')
            xlim([0 12]), ylim([0 12]), xticks([0:6:12]),yticks([0:6:12])
        end

        if iE == 2
            yticklabels('');
        end

        plot(xlim,ylim, ('-k')) % add identity line
        scatter(tmp_subject(:), tmp_model(:), 50, "MarkerFaceColor",plot_colour(iE,:), "MarkerEdgeColor",'w','MarkerFaceAlpha', 0.5,'MarkerEdgeAlpha',0.5)

        [r p] = corr(tmp_subject(:),tmp_model(:), 'type', 'Pearson')

        str=sprintf('r = %1.2f',r);
        T = textbypos(corr_coord(iE),0.87,str);

        if strcmp(study{s}, 'kane')
        ylabel(tl, 'Model SD_{leave} (harvests)')
        xlabel(tl, 'Subject SD_{leave} (harvests)')
        else
        ylabel(tl, 'Model SD_{leave} (s)')
        xlabel(tl, 'Subject SD_{leave} (s)')
        end

    end

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig7_%s_SD_per_participant', study{s})],'-dsvg')
    end

    %% Panel: learning rate fit vs subject LT
    % for winning model

    % get mean leave times across all trials for each subject
    leaveT = readtable(sprintf('%s_trialbytrial.csv', study{s}));
    if strcmp(study{s},'contrerashuerta')
        leaveT = leaveT(leaveT.ben == 1,:); % exclude other condition
    end

    subjMean = zeros(2, numel(unique(leaveT.sub)));
    subjSD = zeros(2, numel(unique(leaveT.sub)));

    for iS = unique(leaveT.sub)'
        for iE = 1:task.nEnviron
            subjMean(iE,iS) = mean(leaveT.leaveT(leaveT.sub == iS & leaveT.env == iE,:));
            subjSD(iE,iS) = std(leaveT.leaveT(leaveT.sub == iS & leaveT.env == iE,:));

        end
    end

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on
    if strcmp(study{s}, 'kane')
    ylabel('Subject harvests per patch');
    else
    ylabel('Subject leaving times (s)');
    end
    xlabel('log(\alpha)');

    scatter(minNLLFitParams.alphaRho, subjMean(1,:), 50, "MarkerFaceColor",color.rich, "MarkerEdgeColor",'w','MarkerFaceAlpha', 0.5,'MarkerEdgeAlpha',0.5), hold on
    scatter(minNLLFitParams.alphaRho, subjMean(2,:),  50, "MarkerFaceColor",color.poor, "MarkerEdgeColor",'w','MarkerFaceAlpha', 0.5,'MarkerEdgeAlpha',0.5)

    set(gca, 'XScale', 'log')

    [r p] = corr([subjMean(1,:),subjMean(2,:)]',[minNLLFitParams.alphaRho; minNLLFitParams.alphaRho], 'type', 'Spearman');
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
        print([export_path, sprintf('fig7_%s_alpha_LT_supplement', study{s})],'-dsvg')
    end

%% Panel: learning rate fit vs subject SD (VARIABILITY)
    % for winning model

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on
    if strcmp(study{s}, 'kane')
    ylabel('Subject harvests per patch');
    else
    ylabel('Subject leaving times (s)');
    end   
    xlabel('log(\alpha)');

    scatter(minNLLFitParams.alphaRho, subjSD(1,:), 50, "MarkerFaceColor",color.rich, "MarkerEdgeColor",'w','MarkerFaceAlpha', 0.5,'MarkerEdgeAlpha',0.5), hold on
    scatter(minNLLFitParams.alphaRho, subjSD(2,:),  50, "MarkerFaceColor",color.poor, "MarkerEdgeColor",'w','MarkerFaceAlpha', 0.5,'MarkerEdgeAlpha',0.5)

    set(gca, 'XScale', 'log')

    [r p] = corr([subjSD(1,:),subjSD(2,:)]',[minNLLFitParams.alphaRho; minNLLFitParams.alphaRho], 'type', 'Spearman');
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
        print([export_path, sprintf('fig7_%s_alpha_SD_supplement', study{s})],'-dsvg')
    end

    %% Panel: beta and rho trajectories - simulations

    load(sprintf('sim_LT_%s_M%d', study{s}, bestModel), '-mat')

    iS = 3; % take an example subject
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on
    plot(trajectories.beta(:,1,iS), 'Color', color.rich, 'LineWidth',1)
    plot(trajectories.beta(:,2,iS), 'Color', color.poor, 'LineWidth',1)
    xlabel('Time in environment (s)'), ylabel('\beta')

    if save_figs == 1
        print([export_path, sprintf('fig7_%s_beta_trajectory', study{s})],'-dsvg')
    end

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on
    plot(trajectories.rho(:,1,iS), 'Color', color.rich, 'LineWidth',1)
    plot(trajectories.rho(:,2,iS), 'Color', color.poor, 'LineWidth',1)
    xlabel('Time in environment (s)'), ylabel('\rho')

    if save_figs == 1
        print([export_path, sprintf('fig7_%s_rho_trajectory', study{s})],'-dsvg')
    end

    clear betas rhos



    % %% Panel: SD in early vs late - simulations
    % 
    % figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    % tl = tiledlayout(1,2, "TileSpacing","tight");
    % 
    % plot_colours = [color.rich; color.poor];
    % corr_coords = [0.25, 0.75];
    % 
    % mean_SD_early = mean(SD_timepoint.early); error_SD_early = std(SD_timepoint.early)./sqrt(size(SD_timepoint.early,1));
    % mean_SD_late = mean(SD_timepoint.late); error_SD_late = std(SD_timepoint.late)./sqrt(size(SD_timepoint.late,1));
    % 
    % subjectMean = mean(subject_leave_times.mean,3); subjectSEM = std(subject_leave_times.mean, [], 3)./sqrt(size(subject_leave_times.mean,3));
    % 
    % for iE = 1:task.nEnviron
    %     % rich environment
    %     nexttile
    %     % connect points
    %     plot([1,2], [SD_timepoint.early(:,iE),SD_timepoint.late(:,iE)], 'Color', [0.7 0.7 0.7]); hold on
    %     % plot points
    %     scatter(1, SD_timepoint.early(:,iE), 70, 'MarkerFaceColor',plot_colours(iE,:), 'MarkerEdgeColor','w')
    %     scatter(2, SD_timepoint.late(:,iE), 70, 'MarkerFaceColor',plot_colours(iE,:), 'MarkerEdgeColor','w')
    %     set(gca,'XLim', [0 3],'XTick', [1:2], 'XTickLabel', {'Early', 'Late'})
    %     errorbar(1, mean_SD_early(iE), error_SD_early(iE), '_', "CapSize", 0,  "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)
    %     errorbar(2, mean_SD_late(iE), error_SD_late(iE), '_', "CapSize", 0,  "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)
    % 
    %     if strcmp(study{s},'kane')
    %         ylim([2.75 4])
    %     elseif strcmp(study{s}, 'contrerashuerta')
    %         ylim([1 10])
    %     elseif strcmp(study{s}, 'leheron')
    %         ylim([2 13])
    %     end
    % 
    % 
    %     [h p] = ttest(SD_timepoint.early(:,iE),SD_timepoint.late(:,iE))
    % 
    %     if p < .001
    %         sig_val = textbypos(corr_coords(iE),0.77,'p < .001');
    %     elseif p >= .001 && p < .01
    %         sig_val = textbypos(corr_coords(iE),0.77,sprintf('p = %1.3f', p));
    %     else
    %         sig_val = textbypos(corr_coords(iE),0.77,sprintf('p = %1.2f', p));
    %     end
    % end
    % 
    % ylabel(tl, 'Model SD_{leave} (s)')
    % tl.Children(1).YTickLabel = [''];
    % 
    % FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    % if save_figs == 1
    %     print([export_path, sprintf('fig7_%s_SD_early_late', study{s})],'-dsvg')
    % end
    % 
    % clear task


    % %% quick plot for presentation
    %         figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on
    %     plot(repelem(0.365,600), 'Color', color.rich, 'LineWidth',2)
    %     plot(repelem(0.38,600), 'Color', color.poor, 'LineWidth',2)
    %     xlabel('Time in environment (s)'), ylabel('\beta')
    %     ylim([0.3 0.45])
    %
    %     if save_figs == 1
    %         print([export_path, sprintf('fig7_%s_beta_trajectory', study{s})],'-dsvg')
    %     end
    %
    %     figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on
    %     plot(repelem(20.6,600), 'Color', color.rich, 'LineWidth',2)
    %     plot(repelem(18.7,600), 'Color', color.poor, 'LineWidth',2)
    %     ylim([16 22])
    %     xlabel('Time in environment (s)'), ylabel('\rho')
    %
    %     if save_figs == 1
    %         print([export_path, sprintf('fig7_%s_rho_trajectory', study{s})],'-dsvg')
    %     end

end
