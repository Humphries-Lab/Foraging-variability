%%% Figure 2 foraging paper - variability
% Emma Scholey 30 Jan 2023
% last updated 11 June 2025

clearvars; close all

addpath(genpath('./plotFunctions'))
addpath('../model/helperFunctions/')
addpath(genpath('../data/experiment_data/'))

run figure_properties_foraging.m

study = {'leheron', 'contrerashuerta', 'kane'};
example_subject = [2 7 5]; % which subjects to plot for examples, for each study

save_figs = 0; % whether to save figures

for s = 1:numel(study)
    task = buildTask(study{s});

    %%% Panel: Coefficient of variation (CV) in same patch and environment
    % 
    % load(sprintf('subject_LT_%s', study{s}), '-mat');
    % 
    % subj_CV = subject_leave_times.sd./subject_leave_times.mean;
    % mean_CV = mean(subj_CV, 3);
    % CV_SEM = std(subj_CV, [], 3)./sqrt(size(subj_CV,3));
    % 
    % figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    % errorbar(1:task.nPatch,mean_CV(1,:),CV_SEM(1,:),lines.exp,'Color',color.rich,'LineWidth',widths.plot); hold on
    % errorbar(1:task.nPatch,mean_CV(2,:),CV_SEM(2,:),lines.exp,'Color',color.poor,'LineWidth',widths.plot)
    % xlabel('Patch yield')
    % ylabel('CV_{leave} (s)')
    % set(gca,'XTick',1:task.nPatch,'XTickLabel',task.patchNames,'XLim',[0 task.nPatch+1],'YLim',[0 0.6])
    % 
    % FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    % if save_figs == 1
    %     print([export_path, sprintf('fig2_%s_subject_CV', study{s})],'-dsvg')
    % end

    %% Panel: SD distributions

    load(sprintf('subject_LT_%s', study{s}), '-mat');

    subj_SD = subject_leave_times.sd;

    plot_colours = [color.rich; color.poor];
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',[20 20 6 10]);
    tl = tiledlayout(3, 2, 'TileSpacing', 'tight');

    ix = 1;
    n_panels = 1:task.nPatch*task.nEnviron;
    for iP = 1:task.nPatch
        for iE = 1:task.nEnviron
           
        ax_handle(iP,iE) = nexttile; hold on

        % x1 = squeeze(subj_SD(1,iP,:)); % rich
        % x2 = squeeze(subj_SD(2,iP,:)); % poor
        % allX   = [x1(:); x2(:)];

        x = squeeze(subj_SD(iE,iP,:));
        % nbins  = 8;                          
        % edges  = linspace(min(x), max(x), nbins+1);
        %histogram(x,edges,'EdgeAlpha', 0.4, 'EdgeColor', 'white', 'FaceColor',plot_colours(iE,:))

        % Kernel density
        [x_dens, f_dens] = ksdensity(x);
        plot(f_dens, x_dens, 'Color', plot_colours(iE,:), 'LineWidth', 2)

        y_fill = [x_dens(1)  x_dens      x_dens(end)];
        x_fill = [0          f_dens(:).' 0];
        fill(x_fill, y_fill, plot_colours(iE,:), ...
            'FaceAlpha',0.3, 'EdgeColor','none');
        
        ix = ix + 1;
        % [x1_dens, f1_dens] = ksdensity(x1);
        % [x2_dens, f2_dens] = ksdensity(x2);
        % 
        % plot(f1_dens, x1_dens, 'Color',color.rich, 'LineWidth',2);
        % plot(f2_dens, x2_dens, 'Color',color.poor, 'LineWidth',2);

        % y1_fill = [x1_dens(1)  x1_dens      x1_dens(end)];
        % x1_fill = [0          f1_dens(:).' 0];
        % fill(x1_fill, y1_fill, color.rich, ...
        %     'FaceAlpha',0.3, 'EdgeColor','none');
        % 
        % y2_fill = [x2_dens(1)  x2_dens      x2_dens(end)];
        % x2_fill = [0          f2_dens(:).' 0];
        % fill(x2_fill, y2_fill, color.poor, ...
        %     'FaceAlpha',0.3, 'EdgeColor','none');
        end
    end

    linkaxes(ax_handle(:))

    ylabel(tl,'Density')
    xlabel(tl,'SD_{leave}')
    %set(gca,'XTick',1:task.nPatch,'XTickLabel',task.patchNames,'XLim',[0 task.nPatch+1],'YLim',[0 0.6])

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    if save_figs == 1
        print([export_path, sprintf('fig2_%s_SD_distribution', study{s})],'-dsvg')
    end

    clear ax_handle
    %% Panel: SD in early vs late

    leaveT = readtable(sprintf('%s_trialbytrial.csv', study{s}));
    if strcmp(study{s},'contrerashuerta')
        leaveT = leaveT(leaveT.ben == 1,:); % exclude other condition
    end
    load(sprintf('%s_subj_var', study{s}), '-mat');

    nSub = numel(unique(leaveT.sub));

    SD_early_env =  zeros([nSub task.nEnviron]);
    SD_late_env  =  zeros([nSub task.nEnviron]);

    SD_early = zeros([nSub 1]);
    SD_late  = zeros([nSub 1]);

    nBlock = task.nBlocks;

    for iS = 1:nSub

        % add block number to leaving times
        subj = leaveT(leaveT.sub == iS,:);
        blockSwitch = subj_var.blockSwitchIndex{iS};
        subj.blockSwitch = blockSwitch(1:length(subj.sub));

        subj.blockNumber = zeros(length(subj.blockSwitch), 1);
        blockNumber = 0;
        for i = 1:length(subj.blockSwitch)
            if subj.blockSwitch(i) == 1
                blockNumber = blockNumber + 1;
            end
            subj.blockNumber(i) = blockNumber;
        end

        % containers for per-block SDs
        early_block_sd     = nan(task.nPatch,1);
        late_block_sd      = nan(task.nPatch,1);
        early_block_sd_env = nan(task.nPatch, task.nEnviron);
        late_block_sd_env  = nan(task.nPatch, task.nEnviron);

        subj.lateTrials = zeros(size(subj.blockNumber));

        for iB = 1:nBlock
            tmp = subj(subj.blockNumber == iB,:);

            % define late trials within this block
            tmp.lateTrials = cumsum(tmp.leaveT) > sum(tmp.leaveT)/2;
            subj.lateTrials(subj.blockNumber == iB) = tmp.lateTrials;
        end

        % within patch SD
        for iP = 1:task.nPatch
            early_trials = subj.leaveT(subj.lateTrials == 0 & subj.patch == iP);
            late_trials  = subj.leaveT(subj.lateTrials == 1 & subj.patch == iP);

            % Compute SD only if n >= 2, otherwise NaN
            if length(early_trials) >= 2
                early_block_sd(iP) = std(early_trials);
            else
                early_block_sd(iP) = NaN;
            end

            if length(late_trials) >= 2
                late_block_sd(iP) = std(late_trials);
            else
                late_block_sd(iP) = NaN;
            end

            %per-environment SDs for this block
            for iE = 1:task.nEnviron
                early_trials_env = subj.leaveT(subj.lateTrials == 0 & subj.patch == iP & subj.env == iE);
                late_trials_env  = subj.leaveT(subj.lateTrials == 1 & subj.patch == iP & subj.env == iE);

                % Compute SD only if n >= 2, otherwise NaN
                if length(early_trials_env) >= 2
                    early_block_sd_env(iP,iE) = std(early_trials_env);
                else
                    early_block_sd_env(iP,iE) = NaN;
                end

                if length(late_trials_env) >= 2
                    late_block_sd_env(iP,iE) = std(late_trials_env);
                else
                    late_block_sd_env(iP,iE) = NaN;
                end
            end

        end

        % average SDs across blocks for this subject
        SD_early(iS) = mean(early_block_sd,'all', 'omitnan');
        SD_late(iS)  = mean(late_block_sd,'all',  'omitnan');

        for iE = 1:task.nEnviron
            SD_early_env(iS,iE) = mean(early_block_sd_env(:,iE), 'omitnan');
            SD_late_env(iS,iE)  = mean(late_block_sd_env(:,iE),  'omitnan');
        end

    end

    mean_SD_early = mean(SD_early); error_SD_early = std(SD_early)./sqrt(nSub);
    mean_SD_late  = mean(SD_late);  error_SD_late  = std(SD_late)./sqrt(nSub);

    mean_SD_early_env = mean(SD_early_env); error_SD_early_env = std(SD_early_env)./sqrt(nSub);
    mean_SD_late_env  = mean(SD_late_env);  error_SD_late_env  = std(SD_late_env)./sqrt(nSub);


    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    tl = tiledlayout(1,2, "TileSpacing","tight");

    plot_colours = [color.rich; color.poor];
    corr_coords = [0.25, 0.75];

    for iE = 1:task.nEnviron
        nexttile
        % connect points
        plot([1,2], [SD_early_env(:,iE),SD_late_env(:,iE)], 'Color', [0.8 0.8 0.8]); hold on
        % plot points
        scatter(1, SD_early_env(:,iE), 50, 'MarkerFaceColor',plot_colours(iE,:), 'MarkerEdgeColor','w', 'LineWidth',1)
        scatter(2, SD_late_env(:,iE), 50, 'MarkerFaceColor',plot_colours(iE,:), 'MarkerEdgeColor','w','LineWidth',1)
        set(gca,'XLim', [0 3],'XTick', [1:2], 'XTickLabel', {'Early', 'Late'})
        errorbar(1, mean_SD_early_env(iE), error_SD_early_env(iE), '_', "CapSize", 0,  "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)
        errorbar(2, mean_SD_late_env(iE), error_SD_late_env(iE), '_', "CapSize", 0,  "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)

        if strcmp(study{s},'kane')
            ylim([0 3])
        elseif strcmp(study{s}, 'contrerashuerta')
            ylim([0 8])
        elseif strcmp(study{s}, 'leheron')
            ylim([0 12])
        end

        %figure, qqplot(CV_late_env(:,iE) - CV_early_env(:,iE))
        [h_sw,p_sw] = swtest(SD_late_env(:,iE) - SD_early_env(:,iE));

        if p_sw < .05 | strcmp(study{s}, 'kane')
            [p, ~, stats] = signrank(SD_early_env(:,iE),SD_late_env(:,iE), 'method','exact')
            meanEffectSize(SD_early_env(:,iE),SD_late_env(:,iE),Paired=true, Effect='robustcohen')

        else
            [~, p, ~, stats] = ttest(SD_early_env(:,iE),SD_late_env(:,iE))
            meanEffectSize(SD_early_env(:,iE),SD_late_env(:,iE),Paired=true, Effect='cohen')

        end

        if p < .001
            sig_val = textbypos(corr_coords(iE),0.77,'p < .001');
        elseif p >= .001 && p < .01
            sig_val = textbypos(corr_coords(iE),0.77,sprintf('p = %1.3f', p));
        else
            sig_val = textbypos(corr_coords(iE),0.77,sprintf('p = %1.2f', p));
        end
    end

    % SD_write = [SD_early, SD_late];
    % writematrix(SD_write, sprintf('../data/experiment_data/%s/%s_SD_early_late.csv', study{s},study{s}))

    if strcmp(study{s},'kane')
        ylabel(tl, 'SD_{leave} (harvests)')
    else
        ylabel(tl, 'SD_{leave} (s)')
    end
    %tl.Children(1).YTickLabel = [''];

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)

    if save_figs == 1
        print([export_path, sprintf('fig2_%s_SD_early_late_env', study{s})],'-dsvg')
    end

    % plotting across environments
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    % connect points
    plot([1,2], [SD_early,SD_late], 'Color', [0.8 0.8 0.8]); hold on
    % plot points
    scatter(1, SD_early, 50, 'MarkerFaceColor',color.scatter, 'MarkerEdgeColor','w', 'LineWidth',1)
    scatter(2, SD_late, 50, 'MarkerFaceColor',color.scatter, 'MarkerEdgeColor','w', 'LineWidth',1)
    set(gca,'XLim', [0 3],'XTick', [1:2], 'XTickLabel', {'Early', 'Late'})
    errorbar(1,mean_SD_early, error_SD_early, '_', "CapSize", 0,  "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)
    errorbar(2, mean_SD_late, error_SD_late, '_', "CapSize", 0,  "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)


    if strcmp(study{s},'kane')
        ylabel('SD_{leave} (harvests)')
    else
        ylabel('SD_{leave} (s)')
    end

    if strcmp(study{s},'kane')
        ylim([0 3])
    elseif strcmp(study{s}, 'contrerashuerta')
        ylim([0 10])
    elseif strcmp(study{s}, 'leheron')
        ylim([0 10])
    end

    %figure, qqplot(CV_late - CV_early)
    [h_sw,p_sw] = swtest(SD_late - SD_early);

    if p_sw < .05 | strcmp(study{s}, 'kane')
        [p, ~, stats] = signrank(SD_early,SD_late, 'method','exact')
        meanEffectSize(SD_early,SD_late,Paired=true, Effect='robustcohen')

        d = SD_early - SD_late;
        medFun = @(data) median(data);
        ci = bootci(5000, medFun, d);   % 5000 bootstrap resamples

    else
        [~, p, ci, stats] = ttest(SD_early,SD_late)
        meanEffectSize(SD_early,SD_late,Paired=true, Effect='cohen')

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
        print([export_path, sprintf('fig2_%s_SD_early_late', study{s})],'-dsvg')
    end


    %% Panel: SD in early vs middle vs late (thirds)

    % SD_early_env =  zeros([nSub task.nEnviron]);
    % SD_mid_env =  zeros([nSub task.nEnviron]);
    % SD_late_env  =  zeros([nSub task.nEnviron]);
    %
    % SD_early = zeros([nSub 1]);
    % SD_mid = zeros([nSub 1]);
    % SD_late  = zeros([nSub 1]);
    %
    % nBlock = task.nBlocks;
    %
    % early_count_trials = zeros([nSub, task.nPatch]);
    % mid_count_trials = zeros([nSub, task.nPatch]);
    % late_count_trials = zeros([nSub, task.nPatch]);
    %
    % for iS = 1:nSub
    %
    %     % add block number to leaving times
    %     subj = leaveT(leaveT.sub == iS,:);
    %     blockSwitch = subj_var.blockSwitchIndex{iS};
    %     subj.blockSwitch = blockSwitch(1:length(subj.sub));
    %
    %     subj.blockNumber = zeros(length(subj.blockSwitch), 1);
    %     blockNumber = 0;
    %     for i = 1:length(subj.blockSwitch)
    %         if subj.blockSwitch(i) == 1
    %             blockNumber = blockNumber + 1;
    %         end
    %         subj.blockNumber(i) = blockNumber;
    %     end
    %
    %     % containers for per-block SDs
    %     early_block_sd     = nan(task.nPatch,1);
    %     mid_block_sd     = nan(task.nPatch,1);
    %     late_block_sd      = nan(task.nPatch,1);
    %     early_block_sd_env = nan(task.nPatch, task.nEnviron);
    %     mid_block_sd_env = nan(task.nPatch, task.nEnviron);
    %     late_block_sd_env  = nan(task.nPatch, task.nEnviron);
    %
    %     subj.lateTrials = zeros(size(subj.blockNumber));
    %
    %     % segment flags: 1 = early, 2 = middle, 3 = late
    %     subj.segment = zeros(size(subj.blockNumber));
    %
    %     for iB = 1:nBlock
    %         tmp = subj(subj.blockNumber == iB,:);
    %
    %         % cumulative sum of leaveT within block
    %         csum  = cumsum(tmp.leaveT);
    %         total = sum(tmp.leaveT);
    %
    %         % thresholds at one-third and two-thirds of total leaveT
    %         thr1 = total/3;
    %         thr2 = 2*total/3;
    %
    %         seg = zeros(size(tmp.leaveT));
    %         seg(csum <= thr1)                    = 1; % early
    %         seg(csum > thr1 & csum <= thr2)      = 2; % middle
    %         seg(csum > thr2)                     = 3; % late
    %
    %         subj.seg(subj.blockNumber == iB) = seg;
    %     end
    %
    %     % within patch SD
    %     for iP = 1:task.nPatch
    %         early_trials = subj.leaveT(subj.seg == 1 & subj.patch == iP);
    %         mid_trials = subj.leaveT(subj.seg == 2 & subj.patch == iP);
    %         late_trials  = subj.leaveT(subj.seg == 3 & subj.patch == iP);
    %
    %         early_count_trials(iS,iP) = length(early_trials);
    %         mid_count_trials(iS,iP) = length(mid_trials);
    %         late_count_trials(iS,iP) = length(late_trials);
    %
    %         % Compute SD only if n >= 2, otherwise NaN
    %         if length(early_trials) >= 2
    %             early_block_sd(iP) = std(early_trials);
    %         else
    %             early_block_sd(iP) = NaN;
    %         end
    %
    %         if length(mid_trials) >= 2
    %             mid_block_sd(iP) = std(mid_trials);
    %         else
    %             mid_block_sd(iP) = NaN;
    %         end
    %
    %         if length(late_trials) >= 2
    %             late_block_sd(iP) = std(late_trials);
    %         else
    %             late_block_sd(iP) = NaN;
    %         end
    %
    %         for iE = 1:task.nEnviron
    %             early_trials_env = subj.leaveT(subj.seg == 1 & subj.patch == iP & subj.env == iE);
    %             mid_trials_env = subj.leaveT(subj.seg == 2 & subj.patch == iP & subj.env == iE);
    %             late_trials_env  = subj.leaveT(subj.seg == 3 & subj.patch == iP & subj.env == iE);
    %
    %             % Compute SD only if n >= 2, otherwise NaN
    %             if length(early_trials_env) >= 2
    %                 early_block_sd_env(iP,iE) = std(early_trials_env);
    %             else
    %                 early_block_sd_env(iP,iE) = NaN;
    %             end
    %
    %             if length(mid_trials_env) >= 2
    %                 mid_block_sd_env(iP,iE) = std(mid_trials_env);
    %             else
    %                 mid_block_sd_env(iP,iE) = NaN;
    %             end
    %
    %             if length(late_trials_env) >= 2
    %                 late_block_sd_env(iP,iE) = std(late_trials_env);
    %             else
    %                 late_block_sd_env(iP,iE) = NaN;
    %             end
    %         end
    %     end
    %
    %     % average SDs across blocks for this subject
    %     SD_early(iS) = mean(early_block_sd,'all', 'omitnan');
    %     SD_mid(iS) = mean(mid_block_sd,'all', 'omitnan');
    %     SD_late(iS)  = mean(late_block_sd,'all',  'omitnan');
    %
    %     for iE = 1:task.nEnviron
    %         SD_early_env(iS,iE) = mean(early_block_sd_env(:,iE), 'omitnan');
    %         SD_mid_env(iS,iE)  = mean(mid_block_sd_env(:,iE),  'omitnan');
    %         SD_late_env(iS,iE)  = mean(late_block_sd_env(:,iE),  'omitnan');
    %     end
    %
    % end
    %
    % mean_SD_early = mean(SD_early); error_SD_early = std(SD_early)./sqrt(nSub);
    % mean_SD_mid = mean(SD_mid); error_SD_mid = std(SD_mid)./sqrt(nSub);
    % mean_SD_late  = mean(SD_late);  error_SD_late  = std(SD_late)./sqrt(nSub);
    %
    % mean_SD_early_env = mean(SD_early_env); error_SD_early_env = std(SD_early_env)./sqrt(nSub);
    % mean_SD_mid_env = mean(SD_early_env); error_SD_mid_env = std(SD_early_env)./sqrt(nSub);
    % mean_SD_late_env  = mean(SD_late_env);  error_SD_late_env  = std(SD_late_env)./sqrt(nSub);
    %
    %
    % figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    % tl = tiledlayout(1,2, "TileSpacing","tight");
    %
    % plot_colours = [color.rich; color.poor];
    %
    % for iE = 1:task.nEnviron
    %     nexttile
    %     % connect points (early -> middle -> late)
    %     plot(1:3, [SD_early_env(:,iE), SD_mid_env(:,iE), SD_late_env(:,iE)], ...
    %         'Color', [0.8 0.8 0.8]); hold on
    %
    %     % plot points
    %     scatter(1, SD_early_env(:,iE), 50, 'MarkerFaceColor',plot_colours(iE,:), ...
    %         'MarkerEdgeColor','w', 'LineWidth',1)
    %     scatter(2, SD_mid_env(:,iE), 50, 'MarkerFaceColor',plot_colours(iE,:), ...
    %         'MarkerEdgeColor','w', 'LineWidth',1)
    %     scatter(3, SD_late_env(:,iE), 50, 'MarkerFaceColor',plot_colours(iE,:), ...
    %         'MarkerEdgeColor','w', 'LineWidth',1)
    %
    %     set(gca,'XLim', [0 4], 'XTick', 1:3, ...
    %         'XTickLabel', {'Early', 'Mid', 'Late'})
    %
    %     % error bars
    %     errorbar(1, mean_SD_early_env(iE), error_SD_early_env(iE), '_', ...
    %         "CapSize", 0, "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)
    %     errorbar(2, mean_SD_mid_env(iE),   error_SD_mid_env(iE), '_', ...
    %         "CapSize", 0, "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)
    %     errorbar(3, mean_SD_late_env(iE),  error_SD_late_env(iE), '_', ...
    %         "CapSize", 0, "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)
    %
    %     if strcmp(study{s},'kane')
    %         ylim([0 3])
    %     elseif strcmp(study{s}, 'contrerashuerta')
    %         ylim([0 8])
    %     elseif strcmp(study{s}, 'leheron')
    %         ylim([0 12])
    %     end
    %
    % end
    %
    % if strcmp(study{s},'kane')
    %     ylabel(tl, 'SD_{leave} (harvests)')
    % else
    %     ylabel(tl, 'SD_{leave} (s)')
    % end
    %
    % FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    % if save_figs == 1
    %     print([export_path, sprintf('fig2_%s_SD_thirds_env', study{s})],'-dsvg')
    % end
    %
    % figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    %
    % % connect points
    % plot(1:3, [SD_early, SD_mid, SD_late]', 'Color', [0.8 0.8 0.8]); hold on
    %
    % % plot points
    % scatter(ones(nSub,1)*1, SD_early, 50, 'MarkerFaceColor',color.scatter, ...
    %     'MarkerEdgeColor','w', 'LineWidth',1)
    % scatter(ones(nSub,1)*2, SD_mid, 50, 'MarkerFaceColor',color.scatter, ...
    %     'MarkerEdgeColor','w', 'LineWidth',1)
    % scatter(ones(nSub,1)*3, SD_late, 50, 'MarkerFaceColor',color.scatter, ...
    %     'MarkerEdgeColor','w', 'LineWidth',1)
    %
    % set(gca,'XLim', [0 4], 'XTick', 1:3, ...
    %     'XTickLabel', {'Early', 'Mid', 'Late'})
    %
    % % error bars on means
    % errorbar(1, mean_SD_early, error_SD_early, '_', ...
    %     "CapSize", 0, "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)
    % errorbar(2, mean_SD_mid,   error_SD_mid,   '_', ...
    %     "CapSize", 0, "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)
    % errorbar(3, mean_SD_late,  error_SD_late,  '_', ...
    %     "CapSize", 0, "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)
    %
    % if strcmp(study{s},'kane')
    %     ylabel('SD_{leave} (harvests)')
    % else
    %     ylabel('SD_{leave} (s)')
    % end
    %
    % if strcmp(study{s},'kane')
    %     ylim([0 3])
    % elseif strcmp(study{s}, 'contrerashuerta')
    %     ylim([0 8])
    % elseif strcmp(study{s}, 'leheron')
    %     ylim([0 12])
    % end
    %
    % FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)

    %% Panel: deviation from optimal for each patch, for example subject (supplementary)

    % leaveT = readtable(sprintf('%s_trialbytrial.csv', study{s}));
    % if strcmp(study{s},'contrerashuerta')
    %     leaveT = leaveT(leaveT.ben == 1,:); % exclude other condition
    % end
    % load(sprintf('%s_subj_var', study{s}), '-mat');
    %
    % figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    % tl = tiledlayout(1,2);
    %
    % iS = example_subject(s);
    % subj = leaveT(leaveT.sub == iS,:);
    % subjBlockSwitch = subj_var.blockSwitchIndex{iS};
    % patchNumber = 1:size(subj,1);
    %
    % % example variability for one subject
    %
    % for iE = 1:task.nEnviron
    %     nexttile
    %     hold on
    %     l = yline(0, '--', 'LineWidth',1, 'Color','k');
    %
    %     tmp = [subj.leaveT(subj.env == iE),subj.patch(subj.env == iE)]; % restrict to single environment
    %
    %     for iP = 1:task.nPatch
    %         lt = tmp(tmp(:,2)==iP,1) - task.optLT(iE,iP);
    %         scatter(find(tmp(:,2)==iP), lt,25, 'MarkerFaceColor', color.patch(iP,:), 'MarkerEdgeColor','w', 'LineWidth',0.5)
    %     end
    %
    %     switch study{s}
    %         case 'leheron'
    %             ylim([-6 6])
    %         case 'contrerashuerta'
    %             ylim([-2 6])
    %         case 'kane'
    %             ylim([-2 8])
    %     end
    %
    %     xlabel('Patch')
    %
    % end
    %
    %
    % if strcmp(study{s}, 'kane')
    %     ylabel(tl, 'Actual - optimal harvests')
    % else
    %     ylabel(tl, 'Actual - optimal leaving time (s)')
    % end
    %
    % FormatFig_For_Export(gcf,fontsize-1,fontname,widths.axis)
    % if save_figs == 1
    %     print([export_path, sprintf('fig2_%s_example', study{s})],'-dsvg')
    % end
    %
    % clear task RR_fig


    % %% Panel: mean LT in early vs late
    % leaveT = readtable(sprintf('%s_trialbytrial.csv', study{s}));
    % if strcmp(study{s},'contrerashuerta')
    %     leaveT = leaveT(leaveT.ben == 1,:); % exclude other condition
    % end
    % load(sprintf('%s_subj_var', study{s}), '-mat');
    % 
    % nSub = numel(unique(leaveT.sub));
    % 
    % LT_early_env =  zeros([nSub task.nEnviron]);
    % LT_late_env  =  zeros([nSub task.nEnviron]);
    % 
    % LT_early = zeros([nSub 1]);
    % LT_late  = zeros([nSub 1]);
    % 
    % nBlock = task.nBlocks;
    % 
    % for iS = 1:nSub
    % 
    %     % add block number to leaving times
    %     subj = leaveT(leaveT.sub == iS,:);
    %     blockSwitch = subj_var.blockSwitchIndex{iS};
    %     subj.blockSwitch = blockSwitch(1:length(subj.sub));
    % 
    %     subj.blockNumber = zeros(length(subj.blockSwitch), 1);
    %     blockNumber = 0;
    %     for i = 1:length(subj.blockSwitch)
    %         if subj.blockSwitch(i) == 1
    %             blockNumber = blockNumber + 1;
    %         end
    %         subj.blockNumber(i) = blockNumber;
    %     end
    % 
    %     % containers for per-block SDs
    %     early_block_lt     = nan(task.nPatch,1);
    %     late_block_lt      = nan(task.nPatch,1);
    %     early_block_lt_env = nan(task.nPatch, task.nEnviron);
    %     late_block_lt_env  = nan(task.nPatch, task.nEnviron);
    % 
    %     subj.lateTrials = zeros(size(subj.blockNumber));
    % 
    %     for iB = 1:nBlock
    %         tmp = subj(subj.blockNumber == iB,:);
    % 
    %         % define late trials within this block
    %         tmp.lateTrials = cumsum(tmp.leaveT) > sum(tmp.leaveT)/2;
    %         subj.lateTrials(subj.blockNumber == iB) = tmp.lateTrials;
    %     end
    % 
    %     % within patch SD
    %     for iP = 1:task.nPatch
    %         early_trials = subj.leaveT(subj.lateTrials == 0 & subj.patch == iP);
    %         late_trials  = subj.leaveT(subj.lateTrials == 1 & subj.patch == iP);
    % 
    %         % Compute SD only if n >= 2, otherwise NaN
    %         if length(early_trials) >= 2
    %             early_block_lt(iP) = mean(early_trials);
    %         else
    %             early_block_lt(iP) = NaN;
    %         end
    % 
    %         if length(late_trials) >= 2
    %             late_block_lt(iP) = mean(late_trials);
    %         else
    %             late_block_lt(iP) = NaN;
    %         end
    % 
    %         %per-environment SDs for this block
    %         for iE = 1:task.nEnviron
    %             early_trials_env = subj.leaveT(subj.lateTrials == 0 & subj.patch == iP & subj.env == iE);
    %             late_trials_env  = subj.leaveT(subj.lateTrials == 1 & subj.patch == iP & subj.env == iE);
    % 
    %             % Compute SD only if n >= 2, otherwise NaN
    %             if length(early_trials_env) >= 2
    %                 early_block_lt_env(iP,iE) = mean(early_trials_env);
    %             else
    %                 early_block_lt_env(iP,iE) = NaN;
    %             end
    % 
    %             if length(late_trials_env) >= 2
    %                 late_block_lt_env(iP,iE) = mean(late_trials_env);
    %             else
    %                 late_block_lt_env(iP,iE) = NaN;
    %             end
    %         end
    % 
    %     end
    % 
    %     % average SDs across blocks for this subject
    %     LT_early(iS) = mean(early_block_lt,'all', 'omitnan');
    %     LT_late(iS)  = mean(late_block_lt,'all',  'omitnan');
    % 
    %     for iE = 1:task.nEnviron
    %         LT_early_env(iS,iE) = mean(early_block_lt_env(:,iE), 'omitnan');
    %         LT_late_env(iS,iE)  = mean(late_block_lt_env(:,iE),  'omitnan');
    %     end
    % 
    % end
    % 
    % mean_early = mean(LT_early); error_early = std(LT_early)./sqrt(nSub);
    % mean_late  = mean(LT_late);  error_late  = std(LT_late)./sqrt(nSub);
    % 
    % mean_early_env = mean(LT_early_env); error_early_env = std(LT_early_env)./sqrt(nSub);
    % mean_late_env  = mean(LT_late_env);  error_late_env  = std(LT_late_env)./sqrt(nSub);
    % 
    % 
    % figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    % tl = tiledlayout(1,2, "TileSpacing","tight");
    % 
    % plot_colours = [color.rich; color.poor];
    % corr_coords = [0.25, 0.75];
    % 
    % for iE = 1:task.nEnviron
    %     ax_handle(iE) = nexttile;
    %     % connect points
    %     plot([1,2], [LT_early_env(:,iE),LT_late_env(:,iE)], 'Color', [0.8 0.8 0.8]); hold on
    %     % plot points
    %     scatter(1, LT_early_env(:,iE), 50, 'MarkerFaceColor',plot_colours(iE,:), 'MarkerEdgeColor','w', 'LineWidth',1)
    %     scatter(2, LT_late_env(:,iE), 50, 'MarkerFaceColor',plot_colours(iE,:), 'MarkerEdgeColor','w','LineWidth',1)
    %     set(gca,'XLim', [0 3],'XTick', [1:2], 'XTickLabel', {'Early', 'Late'})
    %     errorbar(1, mean_early_env(iE), error_early_env(iE), '_', "CapSize", 0,  "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)
    %     errorbar(2, mean_late_env(iE), error_late_env(iE), '_', "CapSize", 0,  "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)
    % 
    %     yline(mean(task.optLT(iE,:)), ':', 'LineWidth', 1.5, 'Color', plot_colours(iE,:))
    % end
    % 
    % linkaxes(ax_handle)
    % 
    % if strcmp(study{s},'kane')
    %     ylabel(tl, 'Leaving time (harvests)')
    % else
    %     ylabel(tl, 'Leaving time (s)')
    % end
    % 
    % if strcmp(study{s},'kane')
    %     ylim([0 10])
    % elseif strcmp(study{s}, 'contrerashuerta')
    %     ylim([0 25])
    % elseif strcmp(study{s}, 'leheron')
    %     ylim([0 50])
    % end
    % 
    % FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    % clear ax_handle
    % % plotting across environments
    % figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    % % connect points
    % plot([1,2], [LT_early,LT_late], 'Color', [0.8 0.8 0.8]); hold on
    % % plot points
    % scatter(1, LT_early, 50, 'MarkerFaceColor',color.scatter, 'MarkerEdgeColor','w', 'LineWidth',1)
    % scatter(2, LT_late, 50, 'MarkerFaceColor',color.scatter, 'MarkerEdgeColor','w', 'LineWidth',1)
    % set(gca,'XLim', [0 3],'XTick', [1:2], 'XTickLabel', {'Early', 'Late'})
    % errorbar(1,mean_early, error_early, '_', "CapSize", 0,  "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)
    % errorbar(2, mean_late, error_late, '_', "CapSize", 0,  "LineWidth", widths.plot, "Color", 'k', 'MarkerSize', 10)
    % 
    % yline(mean(task.optLT,'all'), ':', 'LineWidth', 1.5)
    % 
    % 
    % if strcmp(study{s},'kane')
    %     ylabel('Leaving time (harvests)')
    % else
    %     ylabel('Leaving time (s)')
    % end
    % 
    % if strcmp(study{s},'kane')
    %     ylim([0 10])
    % elseif strcmp(study{s}, 'contrerashuerta')
    %     ylim([0 25])
    % elseif strcmp(study{s}, 'leheron')
    %     ylim([0 50])
    % end
    % 
    % FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)


end
