%%% Figure 3 foraging paper
% Emma Scholey 17 September 2025

clearvars; close all

addpath('./plotFunctions')
addpath(genpath('../data/'))
addpath(genpath('../model/'))

run figure_properties_foraging.m

study = {'leheron', 'contrerashuerta', 'kane'};
save_figs = 1; % whether to save figures

bias_slice = 25; % where to take the slice for the line plot

model = 3; % just plot for one model, since identical
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

% Colourbar across individual patch types
SD_min = min(cellfun(@(x) min(x(:)), SD_all));
SD_max = max(cellfun(@(x) max(x(:)), SD_all));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',[20, 20, 7, 4.5]);
clim([SD_min SD_max]);
colormap(brewermap([], 'Greys'));
colorbar;
axis off
FormatFig_For_Export(gcf,16,fontname,widths.axis)

if save_figs == 1
    print('figures/fig3_SD_colorbar','-dsvg')
end

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
    print('figures/fig3_SD_colorbar_diff','-dsvg')
end


%% DO plotting
for s = 1:numel(study)

    SD_leave = SD_all{s};
    beta_parameter = beta{s};
    bias_parameter = bias{s};

    % -------------- Heatmap - per patch type

    for iR = 1:size(SD_leave, 3) % for each patch
        figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.small_panel);

        % Use imagesc on transformed axes for nice log ticks
        imagesc(log10(beta_parameter), bias_parameter, SD_leave(:,:,iR)); axis xy;

        clim([SD_min SD_max]); %colorbar;
        %colormap(parula);
        colormap(brewermap([], 'Greys'))
        %Beta ticks (log)
        xtick_vals = [1e-3 1e-2 1e-1 1 10];
        xticks(log10(xtick_vals));
        xticklabels(string(xtick_vals));

        hold on

        yline(bias_parameter(bias_slice), 'k--', 'LineWidth',2);

        set(gca,'LineWidth',widths.plot,'Box','off')

        FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)

        if save_figs == 1
            print([sprintf('figures/fig3_SD_heatmap_%s_M%d_patch%d', study{s}, model, iR)],'-dsvg')
        end

    end

    % -------------- Heatmap: difference between high and low yield

    SD_diff = SD_diff_all{s};
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.small_panel);

    imagesc(log10(beta_parameter), bias_parameter, SD_diff); axis xy;

    clim([SD_min_diff SD_max_diff]); %colorbar;
    colormap(brewermap([], 'Greys'))

    %Beta ticks (log)
    xtick_vals = [1e-3 1e-2 1e-1 1 10];
    xticks(log10(xtick_vals));
    xticklabels(string(xtick_vals));

    hold on

    yline(bias_parameter(bias_slice), 'k--', 'LineWidth',2)


    % xlabel('Reward weighting, \beta')
    % ylabel('Bias, c')
    set(gca,'LineWidth',widths.plot,'Box','off')

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)

    if save_figs == 1
        print([sprintf('figures/fig3_SD_heatmap_diff_%s_M%d', study{s}, model)],'-dsvg')
    end

    % ---------- plot slice for one patch and difference

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',[20 20 4.5 3]);



    semilogx(beta_parameter, SD_leave(bias_slice, :, 1), 'LineWidth',widths.plot, 'Color',color.general) % low yield patch
    %ylabel('SD_{leave}'), xlabel('Rew. weight. \beta')

    xticks([0.01 0.1 1]);
    xticklabels({'0.01','0.1','1'});

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)

    if save_figs == 1
        print([sprintf('figures/fig3_SD_heatmap_%s_slice_M%d', study{s}, model)],'-dsvg')
    end

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',[20 20 4.5 3]);

    semilogx(beta_parameter, SD_diff(bias_slice, :, 1), 'LineWidth',widths.plot, 'Color',color.general) % difference
    %ylabel('SD_{leave} diff.'), xlabel('Rew. weight. \beta')
    yline(0,'--', 'LineWidth', 1)

    xticks([0.01 0.1 1]);
    xticklabels({'0.01','0.1','1'});

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)

    if save_figs == 1
        print([sprintf('figures/fig3_SD_heatmap_diff_%s_slice_M%d', study{s}, model)],'-dsvg')
    end

end


