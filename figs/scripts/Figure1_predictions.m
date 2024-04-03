%%% Figure 2 foraging paper - stochastic model expected leaving times
% Emma Scholey 7 Aug 2023

clearvars; close all

addpath('./plotFunctions')
addpath('../../model/helperFunctions/')
addpath('../../data/analytical_data/')

run figure_properties_foraging.m

panelData = {
    'expectedLT_example_linear_rangebias_beta=fixed';...
    'expectedLT_example_linear_rangebeta_bias=-2';...
    'expectedLT_example_exp_rangebias_beta=0';...
    'expectedLT_example_exp_rangebeta_bias=0';...
    'expectedLT_example_exp_rangebias_beta=fixed';...
    'expectedLT_example_exp_rangebeta_bias=-2'
    };

for s = 1:numel(panelData)

    % example reward rate values
    r0=[30 50 70]; % initial yield

    % convert reward rates into optimal leaving times as defined by mvt
    for e=1:2 %each env
        for p=1:3 % each patch type
            if contains(panelData{s}, 'linear')
                optimalAvgRR = [25.84 19.14]; % taken from Nils script with r0 & a as below, and travel time = 5s
                a=4; % decay rate
                optLT(e,p)=(optimalAvgRR(e)-r0(p))/-a;
            else
                optimalAvgRR = [23.40 18.58]; % taken from Nils script with r0 & a as below, and travel time = 5s
                a=0.1; % decay rate
                optLT(e,p)=(log(optimalAvgRR(e)/r0(p)))/-a;
            end
        end
    end


    %% Panel: expected leave times for range of explore parameter 
    load(sprintf('%s.mat', panelData{s}));

    % beta prediction for MVT leaving times
    if contains(panelData{s},'rangebias')
        if strcmp(panelData{s}, 'expectedLT_example_exp_rangebias_beta=0')
            ixRich = find(data.E_leave(:,2) >= optLT(1,2),1,"first");
            ixPoor = find(data.E_leave(:,2) >= optLT(2,2),1,"first");
        else
            ixRich = find(data.E_leave(:,2) >= optLT(1,2),1,"last");
            ixPoor = find(data.E_leave(:,2) >= optLT(2,2),1,"last");
        end
    elseif contains(panelData{s},'rangebeta')
        ixRich = find(data.E_leave(:,2) >= optLT(1,2),1,"first");
        ixPoor = find(data.E_leave(:,2) >= optLT(2,2),1,"first");
    end
    parameter_rich = data.explore(ixRich);
    parameter_poor = data.explore(ixPoor);

    E_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    colororder(color.patch);
    semilogx(data.explore,data.E_leave,'LineWidth',widths.plot); hold on

    if contains(panelData{s},'rangebias')
        xlabel('Bias (higher = exploit)')
        if strcmp(panelData{s}, 'expectedLT_example_exp_rangebias_beta=0')
            xlim([10^-1 10^0.5])
        else
            xlim([-10^2 -10^-1])
        end
    elseif contains(panelData{s},'rangebeta')
        xlabel('Beta (higher = exploit)')
    end

    ylabel('Expected leaving time (s)')
    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    line([parameter_rich parameter_rich],[0 E_fig.Children.YLim(2)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
    line([parameter_poor parameter_poor],[0 E_fig.Children.YLim(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')
    print([export_path 'fig1_' panelData{s}],'-dsvg')
    print([overleaf_path, 'fig1/' 'fig1_' panelData{s}],'-dsvg')

    %% Panel: patch leaving time predictions for model

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
    plot(1:3,optLT(1,:),':','Color',color.rich,'LineWidth',widths.plot); hold on
    plot(1:3,optLT(2,:),':','Color',color.poor,'LineWidth',widths.plot)
    xlabel('Patch yield')
    ylabel('Patch leaving time (s)')
    set(gca,'XTick',1:3,'XTickLabel',{'Low', 'Mid', 'High'},'XLim',[0 4],'YLim',[0 15])

    plot(1:3,data.E_leave(ixRich,:),lines.model,'Color',color.rich,'LineWidth',widths.plot); hold on
    plot(1:3,data.E_leave(ixPoor,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)

    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
    print([export_path 'fig1_patch_times_' panelData{s}],'-dsvg')
    print([overleaf_path, 'fig1/', 'fig1_patch_times_' panelData{s}],'-dsvg')

    clear optLT
end

 %% Panel: exponential reward rates

 r0=[30 50 70]; % initial yield
 optimalAvgRR = [23.40 18.58]; % taken from Nils script with r0 & a, and travel time = 5s, and patch proportions [2/10 3/10 5/10] for rich, [5/10 3/10 2/10] for poor
 a=0.1; % decay rate

 % plot decay function
 RR_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);

 r_t_series = r0.*(exp(-a.*[0:30])');

 colororder(color.patch);
 plot([0:30],r_t_series,'LineWidth', widths.plot)
 %ylim([0 max(30)])

 ylabel('Patch reward rate (units/s)')
 xlabel('Time in patch (s)')
 line([0 RR_fig.Children.XLim(2)],[optimalAvgRR(1) optimalAvgRR(1)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
 line([0 RR_fig.Children.XLim(2)],[optimalAvgRR(2) optimalAvgRR(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')

 FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
 print([export_path, 'fig1_RR_exponential'],'-dsvg')
 print([overleaf_path, 'fig1/', 'fig1_RR_exponential'],'-dsvg')

 %% Panel: linear reward rates

 r0=[30 50 70]; % initial yield
 optimalAvgRR = [25.84 19.14]; % taken from Nils script with r0 & a, and travel time = 5s, and patch proportions [2/10 3/10 5/10] for rich, [5/10 3/10 2/10] for poor
 a=4; % decay rate

 % plot decay function
 RR_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);

 r_t_series = r0-a.*[0:20]';

 colororder(color.patch);
 plot([0:20],r_t_series,'LineWidth', widths.plot)
 ylim([0 80])

 ylabel('Patch reward rate (units/s)')
 xlabel('Time in patch (s)')
 line([0 RR_fig.Children.XLim(2)],[optimalAvgRR(1) optimalAvgRR(1)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
 line([0 RR_fig.Children.XLim(2)],[optimalAvgRR(2) optimalAvgRR(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')

 FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
 print([export_path, 'fig1_RR_linear'],'-dsvg')
 print([overleaf_path, 'fig1/', 'fig1_RR_linear'],'-dsvg')

