%%% Supplementary Fig foraging paper - E leave predictions and subject E
%%% leave
% Emma Scholey 17 September 2025

clearvars; close all

addpath('./plotFunctions')
addpath('../model/helperFunctions/')
addpath(genpath('../data/'))

run figure_properties_foraging.m

study = {'leheron', 'contrerashuerta', 'kane'};
model = [2,3,4];
save_figs = 1; % whether to save figures

%% Panels: simulated vs empirical mean leaving times
for s = 1:numel(study)
    task = buildTask(study{s});
    load(sprintf('subject_LT_%s', study{s}), '-mat');
    subjectMean = mean(subject_leave_times.mean,3); subjectSEM = std(subject_leave_times.mean, [], 3)./sqrt(size(subject_leave_times.mean,3));

    for m = 1:numel(model)

        load(sprintf('fitting_results_M%d_%s',model(m),study{s}), '-mat');
        load(sprintf('sim_LT_%s_M%d', study{s}, model(m)), '-mat')

        % summarise data
        modelMean = mean(simulated_leave_times.mean,3); modelSEM = std(simulated_leave_times.mean, [], 3)./sqrt(size(simulated_leave_times.mean,3));

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
            print([export_path, sprintf('fig_supp_%s_M%d_simulated_mean_LT_overlap', study{s}, model(m))],'-dsvg')
        end


        %% Panel: per-participant fit (simulated vs actual)
        % for winning model, M3

        figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); 
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
            ylabel(tl, 'Model harvests');
        else
            xlabel(tl, 'Subject leaving times (s)');
            ylabel(tl, 'Model leaving times (s)');
        end

        FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
        if save_figs == 1
            print([export_path, sprintf('fig_supp_%s_M%d_per_participant', study{s}, model(m))],'-dsvg')
        end
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
        load(sprintf('expected_LT_%s_%s.mat', study{s}, panelData{p}));
        load(sprintf('subject_LT_%s.mat', study{s}));

        E_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
        colororder(color.patch);

        if strcmp(panelData{p},'rangebias')
            signedlog = @(x) sign(x) .* log10(1 + abs(x));
            data.explore = signedlog(data.explore);
            plot(data.explore, data.E_leave, 'LineWidth', widths.plot); hold on

            xlabel('Bias, c')
            xticks = [-60 -10 -1 0 1 3];  % Choose meaningful bias values
            set(gca, 'XTick', signedlog(xticks))
            set(gca, 'XTickLabel', arrayfun(@num2str, xticks, 'UniformOutput', false))

        elseif strcmp(panelData{p},'rangebeta')
            semilogx(data.explore,data.E_leave,'LineWidth',widths.plot); hold on

            xlabel('Rew. weight, B')
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

        print([export_path 'fig_supp_E_leave_' study{s} '_' panelData{p}],'-dsvg')

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
        print([export_path 'fig_supp_patch_times_' study{s} '_' panelData{p}],'-dsvg')

        clear task.optLT
    end
end