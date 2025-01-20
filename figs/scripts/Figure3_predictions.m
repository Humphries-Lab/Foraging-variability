%%% Figure 2 foraging paper - stochastic model expected leaving times
% Emma Scholey 7 Aug 2023

clearvars; close all

addpath('./plotFunctions')
addpath('../../model/helperFunctions/')
addpath('../../data/analytical_data/')
addpath(genpath('../../data/experiment_data/'))

run figure_properties_foraging.m

study = {'leheron', 'contrerashuerta', 'kane'};

panelData = {
    'rangebeta_bias=-2';...
    'rangebias_beta=0.4';...
    'rangebias_beta=0';...
    };

for s = 1:numel(study)

    task = buildTask(study{s});

    for p = 1:numel(panelData)

        %% Panel: expected leave times for range of exploit parameter
        load(sprintf('expectedLT_%s_%s.mat', study{s}, panelData{p}));
        load(sprintf('subject_LT_%s.mat', study{s}));

        E_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
        colororder(color.patch);
        semilogx(data.explore,data.E_leave,'LineWidth',widths.plot); hold on

        if strcmp(panelData{p},'rangebias_beta=0.4')
            xlabel('{\it c} (higher = exploit)')
            xlim([-10^2 -10^-1])
            if strcmp(study{s}, 'kane')
                xlim([-10^2, -10^0])
            end
        elseif strcmp(panelData{p},'rangebeta_bias=-2')
            xlabel('\beta (higher = exploit)')
            xlim([10^-1.5 10^0])
        elseif strcmp(panelData{p}, 'rangebias_beta=0')
            xlabel('{\it c} (higher = exploit)')
            xlim([10^-0.25 10^0.5])
            xticks([10^-0.2, 10^0.5])
            xticklabels({'10^{-0.2}', '10^{0.5}'})
        end

        if strcmp(study{s}, 'kane')
            ylabel('Expected harvests per patch')
        else
            ylabel('Expected leaving time (s)')
        end


        FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)

        % beta prediction for MVT leaving times, based on mid patch
        if ismember(panelData{p},{'rangebeta_bias=-2', 'rangebias_beta=0'})
            ixRich_MVT = find(data.E_leave(:,2) >= task.optLT(1,2),1,"first");
            ixPoor_MVT = find(data.E_leave(:,2) >= task.optLT(2,2),1,"first");
        elseif ismember(panelData{p}, {'rangebias_beta=0.4'})
            ixRich_MVT = find(data.E_leave(:,2) >= task.optLT(1,2),1,"last");
            ixPoor_MVT = find(data.E_leave(:,2) >= task.optLT(2,2),1,"last");
        end

        parameter_rich = data.explore(ixRich_MVT);
        parameter_poor = data.explore(ixPoor_MVT);
        line([parameter_rich parameter_rich],[0 E_fig.Children.YLim(2)],'Color',[color.rich,0.8], 'LineWidth',widths.plot, 'LineStyle', lines.mvt)
        line([parameter_poor parameter_poor],[0 E_fig.Children.YLim(2)],'Color',[color.poor,0.8], 'LineWidth',widths.plot, 'LineStyle', lines.mvt)

        mean_subj_LT = mean(subject_leave_times.mean,3);
        % beta prediction for subject leaving times, based on mid patch
        if ismember(panelData{p},{'rangebias_beta=0','rangebeta_bias=-2'})
            ixRich_subj = find(data.E_leave(:,2) >= mean_subj_LT(1,2),1,"first");
            ixPoor_subj = find(data.E_leave(:,2) >= mean_subj_LT(2,2),1,"first");
        elseif strcmp(panelData{p}, 'rangebias_beta=0.4')
            ixRich_subj = find(data.E_leave(:,2) >= mean_subj_LT(1,2),1,"last");
            ixPoor_subj = find(data.E_leave(:,2) >= mean_subj_LT(2,2),1,"last");
        end

        parameter_rich = data.explore(ixRich_subj);
        parameter_poor = data.explore(ixPoor_subj);
        line([parameter_rich parameter_rich],[0 E_fig.Children.YLim(2)],'Color',[color.rich,0.8], 'LineWidth',widths.plot)
        line([parameter_poor parameter_poor],[0 E_fig.Children.YLim(2)],'Color',[color.poor,0.8], 'LineWidth',widths.plot)

        print([export_path 'fig3_E_leave_' study{s} '_' panelData{p}],'-dsvg')

        %% Panel: patch leaving time predictions for model

        figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.small_panel);
        % plot(1:3,task.optLT(1,:),lines.mvt,'Color',color.rich,'LineWidth',widths.plot); hold on
        % plot(1:3,task.optLT(2,:),lines.mvt,'Color',color.poor,'LineWidth',widths.plot)

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
