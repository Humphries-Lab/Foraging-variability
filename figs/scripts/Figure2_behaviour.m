%%% Figure 1 foraging paper - behavioural data variability
% Emma Scholey 30 Jan 2023

clearvars; close all

addpath('./plotFunctions')
addpath('../../model/helperFunctions/')
addpath('../../data/experiment_data/')
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

    ylabel('Patch reward rate (units/s)')
    xlabel('Time in patch (s)')
    line([0 RR_fig.Children.XLim(2)],[task.optAvgRR(1) task.optAvgRR(1)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
    line([0 RR_fig.Children.XLim(2)],[task.optAvgRR(2) task.optAvgRR(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig2_%s_RR', study{s})],'-dsvg')
        print([overleaf_path, 'fig2/', sprintf('fig2_%s_RR', study{s})],'-dsvg')
    end

    %% Panel: participant leave times
    load(sprintf('subjectLT_%s', study{s}), '-mat');

    meanLT = mean(subject_leave_times.mean,3);
    subjSEM = std(subject_leave_times.mean, [], 3)./sqrt(size(subject_leave_times.mean,3));

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    errorbar(1:task.nPatch,meanLT(1,:),subjSEM(1,:),lines.exp,'Color',color.rich,'LineWidth',widths.plot); hold on
    errorbar(1:task.nPatch,meanLT(2,:),subjSEM(2,:),lines.exp,'Color',color.poor,'LineWidth',widths.plot)
    xlabel('Patch yield')
    ylabel('Patch leaving time (s)')
    set(gca,'XTick',1:task.nPatch,'XTickLabel',task.patchNames,'XLim',[0 task.nPatch+1],'YLim',[0 round(max(max(meanLT)))+5])

    plot(1:task.nPatch,task.optLT(1,:),':','Color',color.rich,'LineWidth',widths.plot);
    plot(1:task.nPatch,task.optLT(2,:),':','Color',color.poor,'LineWidth',widths.plot)

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig2_%s_subjectLT', study{s})],'-dsvg')
        print([overleaf_path, 'fig2/', sprintf('fig2_%s_subjectLT', study{s})],'-dsvg')
    end

    %% Panel: example participant leaving times across task (supplementary)

    leaveT = readtable(sprintf('%s_trialbytrial.csv', study{s}));
    if strcmp(study{s},'contrerashuerta')
        leaveT = leaveT(leaveT.ben == 1,:); % exclude other condition
    end
    load(sprintf('%s_blockSwitchIndex', study{s}), '-mat');

    optLTenv = mean(task.optLT,2);

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);hold on
    xlabel('Patch number')
    ylabel('Patch leaving time (s)')

    iS = example_subject(s);
    subj = leaveT(leaveT.sub == iS,:);
    subjBlockSwitch = blockSwitchIndex{iS};

    patchNumber = 1:size(subj,1);
    plot(patchNumber, subj.leaveT, 'k')

    blockStarts = find(subjBlockSwitch);
    for iB = 1:task.nBlocks
        if iB == task.nBlocks
            block = patchNumber(blockStarts(iB):end);
        else
            block = patchNumber(blockStarts(iB):blockStarts(iB+1));
        end

        if subj.env(blockStarts(iB)) == 1 % rich
            l = line([block(1), block(end)], [optLTenv(1),optLTenv(1)]); l.Color = color.rich; l.LineWidth = widths.plot;
        elseif subj.env(blockStarts(iB)) == 2 % poor
            l = line([block(1), block(end)], [optLTenv(2),optLTenv(2)]); l.Color = color.poor; l.LineWidth = widths.plot;
        end
    end

    ylim([0 round(max(subj.leaveT))+5])
    legend({'', 'Optimal leave time (rich)', 'Optimal leave time (poor)'}, "Location","southoutside")

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig2_%s_example', study{s})],'-dsvg')
        print([overleaf_path, 'fig2/', sprintf('fig2_%s_example', study{s})],'-dsvg')
    end

    %% Panel: SD in early vs late

    nSub = numel(unique(leaveT.sub));
    earlyTrials = [];
    lateTrials = [];
    for iS = unique(leaveT.sub)'
        subj = leaveT(leaveT.sub == iS,:);
        for iE = 1:task.nEnviron
            tmp = subj.leaveT(subj.env==iE);
            splitIndex = cumsum(tmp) > sum(tmp)/2;
            subj.lateTrials(subj.env == iE) = splitIndex;
        end
        earlyTrials = [earlyTrials; subj(subj.lateTrials == 0,:)];
        lateTrials = [lateTrials; subj(subj.lateTrials == 1,:)];

    end
    total_early_subj_sd =  zeros([nSub 1]);total_late_subj_sd =  zeros([nSub 1]);total_subj_sd =  zeros([nSub 1]);

    for iS = 1:nSub
        total_early_subj_sd(iS) = std(earlyTrials.leaveT(earlyTrials.sub==iS));
        total_late_subj_sd(iS) = std(lateTrials.leaveT(lateTrials.sub==iS));
        total_subj_sd(iS) = std(leaveT.leaveT(leaveT.sub==iS));
    end

    mean(total_subj_sd) % for reporting
    
    early_total_sd= mean(total_early_subj_sd);
    early_total_sd_error = std(total_early_subj_sd)./sqrt(numel(total_early_subj_sd));
    late_total_sd= mean(total_late_subj_sd);
    late_total_sd_error = std(total_late_subj_sd)./sqrt(numel(total_late_subj_sd));
    [h p] = ttest(total_early_subj_sd,total_late_subj_sd); % not significantly different 

    % write to csv for equivalence testing in R
    tmp = [total_early_subj_sd, total_late_subj_sd];
    writematrix(tmp,sprintf('../../data/experiment_data/figures_stats_data/fig2_%s_SD_early_late_data.csv', study{s}))

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    scatter(1, total_early_subj_sd, 70, 'MarkerFaceColor',color.general, 'MarkerEdgeColor','w'); hold on
    scatter(2, total_late_subj_sd, 70, 'MarkerFaceColor',color.general, 'MarkerEdgeColor','w')

    ylabel('SD patch leaving time (s)')
    set(gca,'XLim', [0 3],'XTick', [1:2], 'XTickLabel', {'Early', 'Late'})
    errorbar(1, early_total_sd, early_total_sd_error, '_', "CapSize", 0,  "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)
    errorbar(2, late_total_sd, late_total_sd_error, '_', "CapSize", 0,  "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)
    ylim([round(min([total_early_subj_sd; total_late_subj_sd])) - 1, round(max([total_early_subj_sd; total_late_subj_sd]))+1])

    % connect points
    plot([1,2], [total_early_subj_sd,total_late_subj_sd], 'Color', [0.7 0.7 0.7])

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig2_%s_SD_early_late', study{s})],'-dsvg')
        print([overleaf_path, 'fig2/', sprintf('fig2_%s_SD_early_late', study{s})],'-dsvg')
    end

    clear task RR_fig
end


%% SUPPLEMENT Panel: Leave predictions for range of beta and c
for s = 1:numel(study)
    task = buildTask(study{s});

    % Panel: range beta
    load(sprintf('expectedLT_%s_rangebeta', study{s}), '-mat')

    ixRich = find(data.E_leave(:,2) >= task.optLT(1,2),1,"first");
    ixPoor = find(data.E_leave(:,2) >= task.optLT(2,2),1,"first");

    parameter_rich = data.explore(ixRich);
    parameter_poor = data.explore(ixPoor);

    E_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    colororder(color.patch);
    semilogx(data.explore,data.E_leave,'LineWidth',widths.plot); hold on
    xlabel('Beta (higher = exploit)'), 
    ylabel('Expected leaving time (s)'), ylim([0 40])

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    line([parameter_rich parameter_rich],[0 E_fig.Children.YLim(2)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
    line([parameter_poor parameter_poor],[0 E_fig.Children.YLim(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')
    print([export_path sprintf('fig2_expectedLT_%s_rangebeta_supplement', study{s})],'-dsvg')
    print([overleaf_path sprintf('fig2/fig2_expectedLT_%s_rangebeta_supplement', study{s})],'-dsvg')

    % Panel: range bias
    load(sprintf('expectedLT_%s_rangebias', study{s}), '-mat')

    ixRich = find(data.E_leave(:,2) >= task.optLT(1,2),1,"last");
    ixPoor = find(data.E_leave(:,2) >= task.optLT(2,2),1,"last");

    parameter_rich = data.explore(ixRich);
    parameter_poor = data.explore(ixPoor);

    E_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    colororder(color.patch);
    semilogx(data.explore,data.E_leave,'LineWidth',widths.plot); hold on
    xlabel('Bias (higher = exploit)'), xlim([-10^2 -10^-1])

    ylabel('Expected leaving time (s)')

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    line([parameter_rich parameter_rich],[0 E_fig.Children.YLim(2)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
    line([parameter_poor parameter_poor],[0 E_fig.Children.YLim(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')
    print([export_path sprintf('fig2_expectedLT_%s_rangebias_supplement', study{s})],'-dsvg')
    print([overleaf_path sprintf('fig2/fig2_expectedLT_%s_rangebias_supplement', study{s})],'-dsvg')


end
