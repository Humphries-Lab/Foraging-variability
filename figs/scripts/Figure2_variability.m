%%% Figure 2 foraging paper - variability
% Emma Scholey 30 Jan 2023

clearvars; close all

addpath('./plotFunctions')
addpath('../../model/helperFunctions/')
addpath(genpath('../../data/experiment_data/'))
addpath('../../data/analytical_data/')

run figure_properties_foraging.m

study = {'leheron', 'contrerashuerta', 'kane'};
example_subject = [2 7 5]; % which subjects to plot for examples, for each study

save_figs = 1; % whether to save figures

for s = 1:numel(study)
    task = buildTask(study{s});

    %% Panel: Coefficient of variation (CV) in same patch and environment

    load(sprintf('subject_LT_%s', study{s}), '-mat');

    subj_CV = subject_leave_times.sd./subject_leave_times.mean;
    mean_CV = mean(subj_CV, 3);
    CV_SEM = std(subj_CV, [], 3)./sqrt(size(subj_CV,3));

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    errorbar(1:task.nPatch,mean_CV(1,:),CV_SEM(1,:),lines.exp,'Color',color.rich,'LineWidth',widths.plot); hold on
    errorbar(1:task.nPatch,mean_CV(2,:),CV_SEM(2,:),lines.exp,'Color',color.poor,'LineWidth',widths.plot)
    xlabel('Patch yield')
    ylabel('CV_{leave} (s)')
    set(gca,'XTick',1:task.nPatch,'XTickLabel',task.patchNames,'XLim',[0 task.nPatch+1],'YLim',[0 0.6])

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig2_%s_subject_CV', study{s})],'-dsvg')
    end

    %% Panel: CV in early vs late

    leaveT = readtable(sprintf('%s_trialbytrial.csv', study{s}));
    if strcmp(study{s},'contrerashuerta')
        leaveT = leaveT(leaveT.ben == 1,:); % exclude other condition
    end
    load(sprintf('%s_blockSwitchIndex', study{s}), '-mat');


    nSub = numel(unique(leaveT.sub));

    early_subj_sd_env =  zeros([nSub task.nEnviron]);late_subj_sd_env =  zeros([nSub task.nEnviron]);
    early_subj_mean_env =  zeros([nSub task.nEnviron]);late_subj_mean_env =  zeros([nSub task.nEnviron]);

    early_subj_mean = zeros([nSub 1]); early_subj_sd = zeros([nSub 1]);
    late_subj_mean = zeros([nSub 1]); late_subj_sd = zeros([nSub 1]);

    for iS = 1:nSub

        % add block number to leaving times
        subj = leaveT(leaveT.sub == iS,:);
        blockSwitch = blockSwitchIndex{iS};
        subj.blockSwitch = blockSwitch(1:length(subj.sub));

        subj.blockNumber = zeros(length(subj.blockSwitch), 1);
        blockNumber = 0;
        for i = 1:length(subj.blockSwitch)
            if subj.blockSwitch(i) == 1
                blockNumber = blockNumber + 1;
            end
            subj.blockNumber(i) = blockNumber;
        end

        subj.lateTrials = zeros(size(subj.blockNumber));

        for iB = 1:max(subj.blockNumber)
            tmp = subj(subj.blockNumber == iB,:);
            tmp.lateTrials = cumsum(tmp.leaveT) > sum(tmp.leaveT)/2;

            subj.lateTrials(subj.blockNumber == iB) = tmp.lateTrials;
        end

        early_subj_mean(iS) = mean(subj.leaveT(subj.lateTrials == 0));
        late_subj_mean(iS) = mean(subj.leaveT(subj.lateTrials == 1));
        early_subj_sd(iS) = std(subj.leaveT(subj.lateTrials == 0));
        late_subj_sd(iS) = std(subj.leaveT(subj.lateTrials == 1));

        for iE = 1:task.nEnviron
            early_subj_mean_env(iS,iE) = mean(subj.leaveT(subj.lateTrials == 0 & subj.env == iE));
            late_subj_mean_env(iS,iE) = mean(subj.leaveT(subj.lateTrials == 1 & subj.env == iE));

            early_subj_sd_env(iS,iE) = std(subj.leaveT(subj.lateTrials == 0 & subj.env == iE));
            late_subj_sd_env(iS,iE) = std(subj.leaveT(subj.lateTrials == 1 & subj.env == iE));
        end

    end

    CV_early = early_subj_sd./ early_subj_mean;
    CV_late = late_subj_sd./ late_subj_mean;
    mean_CV_early = mean(CV_early); error_CV_early = std(CV_early)./sqrt(nSub);
    mean_CV_late = mean(CV_late); error_CV_late = std(CV_late)./sqrt(nSub);

    CV_early_env = early_subj_sd_env./early_subj_mean_env;
    CV_late_env = late_subj_sd_env./late_subj_mean_env;
    mean_CV_early_env = mean(CV_early_env); error_CV_early_env = std(CV_early_env)./sqrt(nSub);
    mean_CV_late_env = mean(CV_late_env); error_CV_late_env = std(CV_late_env)./sqrt(nSub);

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    tl = tiledlayout(1,2, "TileSpacing","tight");

    plot_colours = [color.rich; color.poor];
    corr_coords = [0.25, 0.75];

    for iE = 1:task.nEnviron
        nexttile
        % connect points
        plot([1,2], [CV_early_env(:,iE),CV_late_env(:,iE)], 'Color', [0.8 0.8 0.8]); hold on
        % plot points
        scatter(1, CV_early_env(:,iE), 50, 'MarkerFaceColor',plot_colours(iE,:), 'MarkerEdgeColor','w', 'LineWidth',1)
        scatter(2, CV_late_env(:,iE), 50, 'MarkerFaceColor',plot_colours(iE,:), 'MarkerEdgeColor','w','LineWidth',1)
        set(gca,'XLim', [0 3],'XTick', [1:2], 'XTickLabel', {'Early', 'Late'})
        errorbar(1, mean_CV_early_env(iE), error_CV_early_env(iE), '_', "CapSize", 0,  "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)
        errorbar(2, mean_CV_late_env(iE), error_CV_late_env(iE), '_', "CapSize", 0,  "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)

        if strcmp(study{s},'kane')
            ylim([0.3 0.65])
        elseif strcmp(study{s}, 'contrerashuerta')
            ylim([0 0.8])
        elseif strcmp(study{s}, 'leheron')
            ylim([0 0.8])
        end

        %figure, qqplot(CV_late_env(:,iE) - CV_early_env(:,iE))
        [h_sw,p_sw] = swtest(CV_late_env(:,iE) - CV_early_env(:,iE));

        if p_sw < .05 | strcmp(study{s}, 'kane')
            [p, ~, stats] = signrank(CV_early_env(:,iE),CV_late_env(:,iE), 'method','exact')
        else
            [~, p, ~, stats] = ttest(CV_early_env(:,iE),CV_late_env(:,iE))
        end

        if p < .001
            sig_val = textbypos(corr_coords(iE),0.77,'p < .001');
        elseif p >= .001 && p < .01
            sig_val = textbypos(corr_coords(iE),0.77,sprintf('p = %1.3f', p));
        else
            sig_val = textbypos(corr_coords(iE),0.77,sprintf('p = %1.2f', p));
        end
    end

    CV_write = [CV_early, CV_late];   
    %writematrix(CV_write, sprintf('../../data/experiment_data/%s/%s_CV_early_late.csv', study{s},study{s}))

    ylabel(tl, 'CV_{leave} (s)')
    tl.Children(1).YTickLabel = [''];

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig2_%s_CV_early_late_env', study{s})],'-dsvg')
    end

    % plotting across environments
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    % connect points
    plot([1,2], [CV_early,CV_late], 'Color', [0.8 0.8 0.8]); hold on
    % plot points
    scatter(1, CV_early, 50, 'MarkerFaceColor',color.scatter, 'MarkerEdgeColor','w', 'LineWidth',1)
    scatter(2, CV_late, 50, 'MarkerFaceColor',color.scatter, 'MarkerEdgeColor','w', 'LineWidth',1)
    set(gca,'XLim', [0 3],'XTick', [1:2], 'XTickLabel', {'Early', 'Late'})
    errorbar(1,mean_CV_early, error_CV_early, '_', "CapSize", 0,  "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)
    errorbar(2, mean_CV_late, error_CV_late, '_', "CapSize", 0,  "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)


    ylabel('CV_{leave} (s)')

    if strcmp(study{s},'kane')
        ylim([0.4 0.6])
    elseif strcmp(study{s}, 'contrerashuerta')
        ylim([0 0.65])
    elseif strcmp(study{s}, 'leheron')
        ylim([0 0.7])
    end

    %figure, qqplot(CV_late - CV_early)
    [h_sw,p_sw] = swtest(CV_late - CV_early);

    if p_sw < .05 | strcmp(study{s}, 'kane')
        [p, ~, stats] = signrank(CV_early,CV_late, 'method','exact')
    else
        [~, p, ~, stats] = ttest(CV_early,CV_late)
    end

    if p < .001
        sig_val = textbypos(0.5,0.77,'p < .001');
    elseif p >= .001 && p < .01
        sig_val = textbypos(0.5,0.77,sprintf('p = %1.3f', p));
    else
        sig_val = textbypos(0.5,0.77,sprintf('p = %1.2f', p));
    end

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig2_%s_CV_early_late', study{s})],'-dsvg')
    end

%% Panel: deviation from optimal for each patch, for example subject (supplementary)

    leaveT = readtable(sprintf('%s_trialbytrial.csv', study{s}));
    if strcmp(study{s},'contrerashuerta')
        leaveT = leaveT(leaveT.ben == 1,:); % exclude other condition
    end
    load(sprintf('%s_blockSwitchIndex', study{s}), '-mat');

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    tl = tiledlayout(1,2);

    iS = example_subject(s);
    subj = leaveT(leaveT.sub == iS,:);
    subjBlockSwitch = blockSwitchIndex{iS};
    patchNumber = 1:size(subj,1);

    % example variability for one subject

    for iE = 1:task.nEnviron
        nexttile
        hold on
        l = yline(0, '--', 'LineWidth',1, 'Color','k');

        tmp = [subj.leaveT(subj.env == iE),subj.patch(subj.env == iE)]; % restrict to single environment

        for iP = 1:task.nPatch
            lt = tmp(tmp(:,2)==iP,1) - task.optLT(iE,iP);
            scatter(find(tmp(:,2)==iP), lt,25, 'MarkerFaceColor', color.patch(iP,:), 'MarkerEdgeColor','w', 'LineWidth',0.5)
        end

        switch study{s}
            case 'leheron'
                ylim([-6 6])
            case 'contrerashuerta'
                ylim([-2 6])
            case 'kane'
                ylim([-2 7.5])
        end

        xlabel('Patch')

    end


    if strcmp(study{s}, 'kane')
        ylabel(tl, 'Actual - optimal harvests')
    else
        ylabel(tl, 'Actual - optimal leaving time (s)')
    end
    tl.Children(1).YTickLabel = [''];

    FormatFig_For_Export(gcf,fontsize-1,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig2_%s_example', study{s})],'-dsvg')
    end

    clear task RR_fig

end
