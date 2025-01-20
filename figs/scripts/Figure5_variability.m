%%% Figure 4 foraging paper 
% Emma Scholey 7 Aug 2023

clearvars; close all

addpath('./plotFunctions')
addpath('../../data/fitting_data/')
addpath('../../model/helperFunctions/')
addpath(genpath('../../data/experiment_data/'))
addpath('../../data/simulation_data/')
addpath('../../data/analytical_data/')

run figure_properties_foraging.m

study = {'leheron', 'contrerashuerta', 'kane'};

save_figs = 1; % whether to save figures

%% Panel: SD_leave predictions for range of beta 

for s = 1:numel(study)
    task = buildTask(study{s});

    load(sprintf('fitting_results_M3_%s',study{s}), '-mat', 'minNLLFitParams');
    beta_rich = mean(minNLLFitParams.beta_rich);     beta_poor = mean(minNLLFitParams.beta_poor);

    load(sprintf('expectedLT_%s_rangebeta_SD', study{s}), '-mat')

    % Panel: expected SD for range of beta

    SD_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    colororder(color.patch);
    semilogx(data.explore,data.SD_leave,'LineWidth',widths.plot); hold on

    set(gca,"XLim", [10^-3, 10^0])
    xlabel('\beta (higher = exploit)')
    ylabel('SD_{leave}');
    ylim([0 5])
    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    line([beta_rich beta_rich],[0 SD_fig.Children.YLim(2)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
    line([beta_poor beta_poor],[0 SD_fig.Children.YLim(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')

    if save_figs == 1
        print([export_path, sprintf('fig5_%s_expected_SD', study{s})],'-dsvg')
    end

end

%% Panels: simulated vs empirical mean leaving times
% for winning model, M3
for s = 1:numel(study)
    task = buildTask(study{s});

    load(sprintf('fitting_results_M3_%s',study{s}), '-mat');
    load(sprintf('subject_LT_%s', study{s}), '-mat');
    load(sprintf('sim_LT_%s_M3', study{s}), '-mat')

    modelSD = mean(simulated_leave_times.sd,3); modelSEM = std(simulated_leave_times.sd, [], 3)./sqrt(size(simulated_leave_times.sd,3));
    subjectSD = mean(subject_leave_times.sd,3); subjectSEM = std(subject_leave_times.sd, [], 3)./sqrt(size(subject_leave_times.sd,3));

    % plot on one panel to show overlap
     % sd panel
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square)
    % plot the actual experimental data
    errorbar(subjectSD(1,:),subjectSEM(1,:),lines.exp, 'Color',color.rich ,'LineWidth', widths.plot); hold on
    errorbar(subjectSD(2,:),subjectSEM(2,:),lines.exp, 'Color',color.poor,'LineWidth', widths.plot)

    set(gca,'XLim',[0, task.nPatch+1],'XTick',1:task.nPatch,'XTickLabel',task.patchNames)
    
    if strcmp(study{s}, 'kane')
        ylabel('SD_{leave} (harvests)');
    else
        ylabel('SD_{leave} (s)');
    end

    errorbar(modelSD(1,:),modelSEM(1,:),lines.model, 'Color',color.rich ,'LineWidth', widths.plot)
    errorbar(modelSD(2,:),modelSEM(2,:),lines.model, 'Color',color.poor,'LineWidth', widths.plot)

    if strcmp(study{s},'kane') % set manually for Kane since subjects SD so similar
         ylim([1 3])
    elseif strcmp(study{s}, 'contrerashuerta')
        ylim([1 4])
    elseif strcmp(study{s}, 'leheron')
        ylim([2 6])
    end

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig5_%s_simulated_sd_LT_overlap', study{s})],'-dsvg')
    end

    % calculate coefficient of variation
    subjCV = mean(subject_leave_times.sd ./ subject_leave_times.mean, 'all');

    %% Panel: per-participant fit (simulated vs actual)
    % for winning model, M3

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on
    tl = tiledlayout(1,2,'TileSpacing','tight');

    plot_colour = [color.rich;color.poor];

    corr_coord = [0.27,0.63];

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
            corr_coord (1) = 0.3; 
        end

        if iE == 2
           yticklabels('');
        end
        
        plot(xlim,ylim, ('-k')) % add identity line
        scatter(tmp_subject(:), tmp_model(:), 50, "MarkerFaceColor",plot_colour(iE,:), "MarkerEdgeColor",'w','MarkerFaceAlpha', 0.5,'MarkerEdgeAlpha',0.5)

        if strcmp(study{s}, 'kane')
        ylabel(tl, 'Model SD_{leave} (harvests)')
        xlabel(tl, 'Subject SD_{leave} (harvests)')
        else
        ylabel(tl, 'Model SD_{leave} (s)')
        xlabel(tl, 'Subject SD_{leave} (s)')
        end

        [r p] = corr(tmp_subject(:),tmp_model(:), 'type', 'Pearson')

        str=sprintf('r = %1.2f',r);
        T = textbypos(corr_coord(iE),0.87,str);

        if p < .001
            sig_val = textbypos(corr_coord(iE),0.77,'p < .001');
        elseif p >= .001 && p < .01
            sig_val = textbypos(corr_coord(iE),0.77,sprintf('p = %1.3f', p));
        else
            sig_val = textbypos(corr_coord(iE),0.77,sprintf('p = %1.2f', p));
        end

    end

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig5_%s_SD_per_participant', study{s})],'-dsvg')
    end

end

