%%% Figure 4 foraging paper
% Emma Scholey 17 September 2025

clearvars; close all

addpath('./plotFunctions')
addpath(genpath('../data/'))
addpath(genpath('../model/'))

run figure_properties_foraging.m

study = {'leheron', 'contrerashuerta', 'kane'};
save_figs = 1; % whether to save figures

model = 4; 
model_table = readtable("foragingModelTable.xlsx");
model_params = model_table(model_table.modelNumber == model,:);

SD_all = cell(numel(study), 1);
SD_diff_all = cell(numel(study), 1);
beta = cell(numel(study), 1);
bias = cell(numel(study), 1);

% load all studies to get the overall max/min SD (align colourbar)
for s = 1:numel(study)
    load(sprintf("expected_SD_%s.mat", study{s}));

    SD_all{s} = data.SD_leave_heatmap;
    SD_diff_all{s} = data.SD_leave_heatmap_diff;
    beta{s} = data.beta;
    bias{s} = data.bias;
end

%% Create stand-alone colorbar across the studies
% Colourbar for differences in SD_leave
SD_min_diff = min(cellfun(@(x) min(x(:)), SD_diff_all));
SD_max_diff = max(cellfun(@(x) max(x(:)), SD_diff_all));

% Create stand-alone colorbar across the studies
figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',[20, 20, 7, 4.5]);
clim([SD_min_diff SD_max_diff]);
colormap(brewermap([], 'Greys'));
colorbar;
axis off
FormatFig_For_Export(gcf,16,fontname,widths.axis)

if save_figs == 1
    print('figures/fig4_SD_colorbar_diff','-dsvg')
end


for s = 1:numel(study)

    beta_parameter = beta{s};
    bias_parameter = bias{s};
    SD_diff = SD_diff_all{s};
    SD_leave = SD_all{s};


    % -------------- Load participants model fits
    load(sprintf('../data/fitting_data/fitting_results_M%d_%s.mat', model, study{s}), 'minNLLFitParams');

    if strcmp(model_params.betaFunction, 'separate')
        beta_rich_fit = minNLLFitParams.beta_rich;
        beta_poor_fit = minNLLFitParams.beta_poor;

        % Precompute summaries
        beta_rich_prc = prctile(beta_rich_fit, [25 50 75]);
        beta_poor_prc = prctile(beta_poor_fit, [25 50 75]);

        % Transform beta to log10 to match your imagesc x-axis
        beta_rich_fit_log = log10(beta_rich_fit);
        beta_rich_prc_log = log10(beta_rich_prc);
        beta_poor_fit_log = log10(beta_poor_fit);
        beta_poor_prc_log = log10(beta_poor_prc);

        % IQR box corners
        x0_r = beta_rich_prc_log(1);  x1_r = beta_rich_prc_log(3);
        x0_p = beta_poor_prc_log(1);  x1_p = beta_poor_prc_log(3);

        % Coordinates of the rectangle corners
        X_r = [x0_r x1_r x1_r x0_r];
        X_p = [x0_p x1_p x1_p x0_p];

    elseif strcmp(model_params.betaFunction, 'single')
        beta_fit = minNLLFitParams.beta;

        beta_prc = prctile(beta_fit, [25 50 75]);

        beta_prc_log = log10(beta_prc);
        beta_fit_log = log10(beta_fit);

        x0 = beta_prc_log(1);  x1 = beta_prc_log(3);

        X = [x0 x1 x1 x0];

    end

    if strcmp(model_params.bias, 'separate')
        bias_rich_fit = minNLLFitParams.bias_rich;
        bias_poor_fit = minNLLFitParams.bias_poor;

        % Precompute summaries
        bias_rich_prc = prctile(bias_rich_fit, [25 50 75]);
        bias_poor_prc = prctile(bias_poor_fit, [25 50 75]);

        % IQR box corners
        y0_r = bias_rich_prc(1);  y1_r = bias_rich_prc(3);
        y0_p = bias_poor_prc(1);  y1_p = bias_poor_prc(3);

        % Coordinates of the rectangle corners
        Y_r = [y0_r y0_r y1_r y1_r];
        Y_p = [y0_p y0_p y1_p y1_p];

    elseif strcmp(model_params.bias, 'single')
        bias_fit = minNLLFitParams.bias;

        bias_prc = prctile(bias_fit, [25 50 75]);

        y0 = bias_prc(1);  y1 = bias_prc(3);

        Y = [y0 y0 y1 y1];

    end


    % -------------- Heatmap: difference between high and low yield

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.small_panel);

    imagesc(log10(beta_parameter), bias_parameter, SD_diff); axis xy;

    clim([SD_min_diff SD_max_diff]); %colorbar;
    colormap(brewermap([], 'Greys'))

    %Beta ticks (log)
    xtick_vals = [1e-3 1e-2 1e-1 1 10];
    xticks(log10(xtick_vals));
    xticklabels(string(xtick_vals));

    hold on
    if strcmp(model_params.betaFunction,'separate') && strcmp(model_params.bias,'single')
        patch(X_r, Y, color.rich, 'FaceAlpha',0,'EdgeColor',color.rich,'LineWidth',0.7);
        patch(X_p, Y, color.poor, 'FaceAlpha',0,'EdgeColor',color.poor,'LineWidth',0.7);
    elseif strcmp(model_params.betaFunction,'single') && strcmp(model_params.bias,'separate')
        patch(X, Y_r, color.rich, 'FaceAlpha',0,'EdgeColor',color.rich,'LineWidth',0.7);
        patch(X, Y_p, color.poor, 'FaceAlpha',0,'EdgeColor',color.poor,'LineWidth',0.7);
    elseif strcmp(model_params.betaFunction,'separate') && strcmp(model_params.bias,'separate')
        patch(X_r, Y_r, color.rich, 'FaceAlpha',0,'EdgeColor',color.rich,'LineWidth',0.7);
        patch(X_p, Y_p, color.poor, 'FaceAlpha',0,'EdgeColor',color.poor,'LineWidth',0.7);
    end

    % xlabel('Reward weighting, \beta')
    % ylabel('Bias, c')
    set(gca,'LineWidth',widths.plot,'Box','off')

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)

    if save_figs == 1
        print([sprintf('figures/fig4_SD_heatmap_diff_%s_M%d_params', study{s}, model)],'-dsvg')
    end

    % ---------- plot slice for all patches

    if strcmp(model_params.betaFunction,'separate') && strcmp(model_params.bias,'single')
        figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',[20 20 4.5 4]);

        median_bias = median(minNLLFitParams.bias);
        median_beta_rich = median(minNLLFitParams.beta_rich);
        median_beta_poor = median(minNLLFitParams.beta_poor);

        bias_ix = find(bias_parameter >= median_bias, 1, 'first');

        colororder(color.patch)
        semilogx(beta_parameter, squeeze(SD_leave(bias_ix, :, :)), 'LineWidth',widths.plot)
        xline(median_beta_rich, 'LineWidth', widths.plot, 'Color',color.rich, 'LineStyle','--')
        xline(median_beta_poor, 'LineWidth', widths.plot, 'Color',color.poor, 'LineStyle','--')

        %xline(median_beta_rich, 'LineWidth', widths.plot, 'LineStyle', ':')
        %xline(median_beta_poor, 'LineWidth', widths.plot)

        xticks([0.01 0.1 1]);
        xticklabels({'0.01','0.1','1'});


        FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
        if save_figs == 1
            print([sprintf('figures/fig4_SD_heatmap_%s_slice_M%d_params_patch_beta', study{s}, model)],'-dsvg')
        end

    elseif strcmp(model_params.betaFunction,'single') && strcmp(model_params.bias,'separate')
        figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',[20 20 2 4.5]);

        median_beta = median(minNLLFitParams.beta);
        median_bias_rich = median(minNLLFitParams.bias_rich);
        median_bias_poor = median(minNLLFitParams.bias_poor);

        beta_ix = find(beta_parameter >= median_beta, 1, 'first');

        colororder(color.patch)
        plot(squeeze(SD_leave(:, beta_ix, :)), bias_parameter, 'LineWidth',widths.plot)

        yline(median_bias_rich, 'LineWidth', widths.plot, 'Color',color.rich, 'LineStyle','--')
        yline(median_bias_poor, 'LineWidth', widths.plot, 'Color',color.poor, 'LineStyle','--')
       % yline(median_bias_rich, 'LineWidth', widths.plot,'LineStyle',':')
       % yline(median_bias_poor, 'LineWidth', widths.plot)

        ylim([-inf, 1])
        xlim([0, inf])
        increment = round(max(max(squeeze(SD_leave(:, beta_ix, :))))/2);
        xticks(0:increment:8)

        FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
        if save_figs == 1
            print([sprintf('figures/fig4_SD_heatmap_%s_slice_M%d_params_patch_bias', study{s}, model)],'-dsvg')
        end

    elseif strcmp(model_params.betaFunction,'separate') && strcmp(model_params.bias,'separate')

        figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',[20 20 4.5 5]);
        t = tiledlayout(2,1, "TileSpacing","compact");

        median_bias_rich = median(minNLLFitParams.bias_rich);
        median_bias_poor = median(minNLLFitParams.bias_poor);

        median_beta_rich = median(minNLLFitParams.beta_rich);
        median_beta_poor = median(minNLLFitParams.beta_poor);

        bias_ix_rich = find(bias_parameter >= median_bias_rich, 1, 'first');
        bias_ix_poor = find(bias_parameter >= median_bias_poor, 1, 'first');
        beta_ix_rich = find(beta_parameter >= median_beta_rich, 1, 'first');
        beta_ix_poor = find(beta_parameter >= median_beta_poor, 1, 'first');

        colororder(color.patch)
        ax1 = nexttile;
        semilogx(beta_parameter, squeeze(SD_leave(bias_ix_rich, :, :)), 'LineWidth',widths.plot)

        xline(median_beta_rich, 'LineWidth', widths.plot, 'Color',color.rich, 'LineStyle','--')
        xline(median_beta_poor, 'LineWidth', widths.plot, 'Color',color.poor, 'LineStyle','--')

        %xline(median_beta_rich, 'LineWidth', widths.plot, 'LineStyle', ':')
        %xline(median_beta_poor, 'LineWidth', widths.plot)


         ax2 =  nexttile;
        semilogx(beta_parameter, squeeze(SD_leave(bias_ix_poor, :, :)), 'LineWidth',widths.plot)
        
        xline(median_beta_rich, 'LineWidth', widths.plot, 'Color',color.rich, 'LineStyle','--')
        xline(median_beta_poor, 'LineWidth', widths.plot, 'Color',color.poor, 'LineStyle','--')

        %xline(median_beta_rich, 'LineWidth', widths.plot, 'LineStyle', ':')
        %xline(median_beta_poor, 'LineWidth', widths.plot)

        xticks([0.01 0.1 1]);
        xticklabels({'0.01','0.1','1'});

        linkaxes([ax1, ax2], 'x');  % Share x-axis
        % Remove x-axis tick labels from the top plot
        ax1.XTickLabel = [];

        FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
        if save_figs == 1
            print([sprintf('figures/fig4_SD_heatmap_%s_slice_M%d_params_patch_beta', study{s}, model)],'-dsvg')
        end

        figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',[20 20 4 4.5]);
        t = tiledlayout(1,2, 'TileSpacing','compact');

        colororder(color.patch)
        ax1 = nexttile;
        plot(squeeze(SD_leave(:, beta_ix_rich, :)), bias_parameter, 'LineWidth',widths.plot)

        %yline(median_bias_rich, 'LineWidth', widths.plot,'LineStyle',':')
        %yline(median_bias_poor, 'LineWidth', widths.plot)

        yline(median_bias_rich, 'LineWidth', widths.plot, 'Color',color.rich, 'LineStyle','--')
        yline(median_bias_poor, 'LineWidth', widths.plot, 'Color',color.poor, 'LineStyle','--')

        ylim([-inf, 1])
        xlim([0, inf])

        ax2 = nexttile;
        plot(squeeze(SD_leave(:, beta_ix_poor, :)), bias_parameter,  'LineWidth',widths.plot)

        %yline(median_bias_rich, 'LineWidth', widths.plot,'LineStyle',':')
        %yline(median_bias_poor, 'LineWidth', widths.plot)
        yline(median_bias_rich, 'LineWidth', widths.plot, 'Color',color.rich, 'LineStyle','--')
        yline(median_bias_poor, 'LineWidth', widths.plot, 'Color',color.poor, 'LineStyle','--')

        ylim([-inf, 1])
        xlim([0, inf])

        linkaxes([ax1 ax2],'xy')
        ax2.YTickLabel = [];

        FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
        if save_figs == 1
            print([sprintf('figures/fig4_SD_heatmap_%s_slice_M%d_params_patch_bias', study{s}, model)],'-dsvg')
        end

    end

%% Panels: simulated vs empirical SD leaving times
    task = buildTask(study{s});
    load(sprintf('subject_LT_%s', study{s}), '-mat');
    subjectSD = mean(subject_leave_times.sd,3); subjectSEM = std(subject_leave_times.sd, [], 3)./sqrt(size(subject_leave_times.sd,3));

    for m = 1:numel(model)

        load(sprintf('fitting_results_M%d_%s',model(m),study{s}), '-mat');
        load(sprintf('sim_LT_%s_M%d', study{s}, model(m)), '-mat')

        modelSD = mean(simulated_leave_times.sd,3); modelSEM = std(simulated_leave_times.sd, [], 3)./sqrt(size(simulated_leave_times.sd,3));

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
            ylim([1 2.5])
        elseif strcmp(study{s}, 'contrerashuerta')
            ylim([1 4])
        elseif strcmp(study{s}, 'leheron')
            ylim([2 6])
        end

        FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
        if save_figs == 1
            print([export_path, sprintf('fig4_%s_M%d_simulated_sd_LT_overlap', study{s}, model(m))],'-dsvg')
        end

        % calculate coefficient of variation
        subjCV = mean(subject_leave_times.sd ./ subject_leave_times.mean, 'all');

        %% Panel: per-participant fit SD (simulated vs actual)

        figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); 
        tl = tiledlayout(1,2,'TileSpacing','tight');

        plot_colour = [color.rich;color.poor];

        corr_coord = [0.27,0.63];

        for iE = 1:task.nEnviron
            nexttile, hold on

            tmp_subject = squeeze(subject_leave_times.sd(iE,:,:));
            tmp_model = squeeze(simulated_leave_times.sd(iE,:,:));

            if strcmp(study{s},'kane')
                xlim([1 2.5]), ylim([1 2.5]), xticks([1:1:3]),yticks([1:1:3])
            elseif strcmp(study{s}, 'contrerashuerta')
                xlim([0 7]), ylim([0 7]), xticks([0:4:8]),yticks([0:4:8])
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
                ylabel(tl, 'Model SD_{leave}')
                xlabel(tl, 'Subject SD_{leave} (harvests)')
            else
                ylabel(tl, 'Model SD_{leave} (s)')
                xlabel(tl, 'Subject SD_{leave} (s)')
            end

            [r p] = corr(tmp_subject(:),tmp_model(:), 'type', 'Pearson');

            str=sprintf('r = %1.2f',r);
            T = textbypos(corr_coord(iE),0.87,str);

            % if p < .001
            %     sig_val = textbypos(corr_coord(iE),0.77,'p < .001');
            % elseif p >= .001 && p < .01
            %     sig_val = textbypos(corr_coord(iE),0.77,sprintf('p = %1.3f', p));
            % else
            %     sig_val = textbypos(corr_coord(iE),0.77,sprintf('p = %1.2f', p));
            % end

        end

        FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
        if save_figs == 1
            print([export_path, sprintf('fig4_%s_M%d_SD_per_participant', study{s}, model(m))],'-dsvg')
        end

    end
end


