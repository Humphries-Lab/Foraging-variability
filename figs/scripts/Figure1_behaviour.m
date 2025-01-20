%%% Figure 1 foraging paper - behavioural and task data
% Emma Scholey 30 Jan 2023

clearvars; close all

addpath('./plotFunctions')
addpath('../../model/helperFunctions/')
addpath(genpath('../../data/experiment_data/'))
addpath('../../data/analytical_data/')

run figure_properties_foraging.m

study = {'leheron', 'contrerashuerta', 'kane'};
example_subject = [2 7 4]; % which subjects to plot for examples, for each study

save_figs = 1; % whether to save figures

for s = 1:numel(study)
    task = buildTask(study{s});

    %% Panel: reward rates

    % plot decay function
    RR_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);

    if strcmp(task.rewardFunction, 'exponential')
        r_t_series = task.r0.*(exp(-task.decayRate.*[0:30])');
    elseif strcmp(task.rewardFunction, 'linear')
        r_t_series = task.r0-task.decayRate.*[0:20]';
    end

    colororder(color.patch);
    if strcmp(task.rewardFunction, 'exponential')
        plot([0:30],r_t_series,'LineWidth', widths.plot)
        ylim([0 max(task.r0+10)])

    elseif strcmp(task.rewardFunction, 'linear')
        plot([0:20],r_t_series,'LineWidth', widths.plot)
        ylim([0 max(task.r0+20)])

    end

    ylabel('Patch reward rate (r/s)')
    if strcmp(study{s}, 'kane')
        xlabel('Harvests in patch')
    else
        xlabel('Time in patch (s)')
    end

    line([0 RR_fig.Children.XLim(2)],[task.optAvgRR(1) task.optAvgRR(1)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
    line([0 RR_fig.Children.XLim(2)],[task.optAvgRR(2) task.optAvgRR(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig1_%s_RR', study{s})],'-dsvg')
    end

    %% Panel: participant leave times
    load(sprintf('subject_LT_%s', study{s}), '-mat');

    meanLT = mean(subject_leave_times.mean,3);
    subjSEM = std(subject_leave_times.mean, [], 3)./sqrt(size(subject_leave_times.mean,3));

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    errorbar(1:task.nPatch,meanLT(1,:),subjSEM(1,:),lines.exp,'Color',color.rich,'LineWidth',widths.plot); hold on
    errorbar(1:task.nPatch,meanLT(2,:),subjSEM(2,:),lines.exp,'Color',color.poor,'LineWidth',widths.plot)
    xlabel('Patch yield')

    if strcmp(study{s}, 'kane')
        ylabel('Harvests per patch')
    else
        ylabel('Patch leaving time (s)')
    end

    set(gca,'XTick',1:task.nPatch,'XTickLabel',task.patchNames,'XLim',[0 task.nPatch+1],'YLim',[0 round(max(max(meanLT)))+5])

    plot(1:task.nPatch,task.optLT(1,:),lines.mvt,'Color',color.rich,'LineWidth',widths.plot);
    plot(1:task.nPatch,task.optLT(2,:),lines.mvt,'Color',color.poor,'LineWidth',widths.plot)

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig1_%s_subject_LT', study{s})],'-dsvg')
    end

end
