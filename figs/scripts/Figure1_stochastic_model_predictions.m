%%% Figure 1 foraging paper - stochastic model expected leaving times 
% Emma Scholey 7 Aug 2023 

clearvars; close all

dataPath = '../../data/fig_data/fig1/';

addpath('./plotFunctions')

run figure_properties_foraging.m

%% Panel: deterministic MVT choice - avgRR vs patch reward rate (exponential decay function)

% example reward rate values (taken from Le Heron) 
optimalAvgRR = [21.8678 18.5632]; 
r0=[32.5 45 57.5]; % initial yield 
a=0.075; % decay rate

 % plot decay function
RR_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);

r_t_series = r0.*(exp(-a.*[1:30])');

colororder(color.patch);
plot(r_t_series,'LineWidth', widths.plot)
ylabel('Patch reward rate (units/s)')
xlabel('Time in patch (s)')

line([0 RR_fig.Children.XLim(2)],[optimalAvgRR(1) optimalAvgRR(1)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', ':')
line([0 RR_fig.Children.XLim(2)],[optimalAvgRR(2) optimalAvgRR(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', ':')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig1_MVT_RR_deterministic'],'-dsvg')


%% Panel: optimal leaving times for MVT 

% convert reward rates into optimal leaving times as defined by mvt
for e=1:2 %each env
    for p=1:3 % each patch type
        optLT(e,p)=(log(optimalAvgRR(e)/r0(p)))/-a; 
    end
end

MVT_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
plot(1:3,optLT(1,:),lines.exp,'Color',color.rich,'LineWidth',1.5); hold on
plot(1:3,optLT(2,:),lines.exp,'Color',color.poor,'LineWidth',1.5)
xlabel('Patch yield')
ylabel('Patch leaving time (s)')
set(gca,'XTick',1:3,'XTickLabel',{'Low', 'Medium', 'High'},'XLim',[0 4],'YLim',[0 30])

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig1_MVT_LT'],'-dsvg')

%% Panel: expected leave times for range of beta (softmax model)

load([dataPath 'expectedLT_softmax_leheron'])

% beta prediction for MVT
ixRich = find(round(data.E_leave(:,2),1) >= round(optLT(1,2),1),1,"first");
ixPoor = find(round(data.E_leave(:,2),1) >= round(optLT(2,2),1),1,"first");

beta_rich = data.beta(ixRich);
beta_poor = data.beta(ixPoor);

E_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
colororder(color.patch);
semilogx(data.beta,data.E_leave,'LineWidth',widths.plot); hold on
xlabel('Beta (higher = exploit)')
ylabel('Expected leaving time (s)')
line([beta_rich beta_rich],[0 E_fig.Children.YLim(2)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
line([beta_poor beta_poor],[0 E_fig.Children.YLim(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig1_LT_beta_range_softmax'],'-dsvg')

%% Panel: expected leave times vs. optimal MVT 

figure(MVT_fig)
plot(1:3,data.E_leave(ixRich,:),lines.model,'Color',color.rich,'LineWidth',widths.plot); hold on
plot(1:3,data.E_leave(ixPoor,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig1_MVT_expectedLT_softmax'],'-dsvg')

%% Panel: expected leave times for range of beta when beta = overharvesting (softmax model) 

beta_rich = 0.20; ixRich = find(round(data.beta,2) == beta_rich,1,"first");
beta_poor = 0.25; ixPoor = find(round(data.beta,2) == beta_poor,1,"first");

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
colororder(color.patch);
semilogx(data.beta,data.E_leave,'LineWidth',widths.plot); hold on
xlabel('Beta (higher = exploit)')
ylabel('Expected leaving time (s)')

% plot overharvesting betas 
line([beta_rich beta_rich],[0 E_fig.Children.YLim(2)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
line([beta_poor beta_poor],[0 E_fig.Children.YLim(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig1_LT_beta_range_softmax_overharvesting'],'-dsvg')

%% Panel: expected overharvesting leave times vs. optimal MVT 

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
plot(1:3,optLT(1,:),lines.exp,'Color',color.rich,'LineWidth',widths.plot); hold on
plot(1:3,optLT(2,:),lines.exp,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('Patch leaving time (s)')
set(gca,'XTick',1:3,'XTickLabel',{'Low', 'Medium', 'High'},'XLim',[0 4],'YLim',[0 30])

plot(1:3,data.E_leave(ixRich,:),lines.model,'Color',color.rich,'LineWidth',widths.plot); hold on
plot(1:3,data.E_leave(ixPoor,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig1_MVT_expectedLT_softmax_overharvesting'],'-dsvg')

%% Panel: expected standard deviation of leave times for range of beta when beta = overharvesting (softmax model) 

beta_rich = 0.22; ixRich = find(round(data.beta,2) == beta_rich,1,"first");
beta_poor = 0.28; ixPoor = find(round(data.beta,2) == beta_poor,1,"first");

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
colororder(color.patch);
semilogx(data.beta,data.SD_leave,'LineWidth',widths.plot); hold on
xlabel('Beta (higher = exploit)')
ylabel('SD leaving times (s)')
ylim([0 10])
% plot overharvesting betas 
line([beta_rich beta_rich],[0 E_fig.Children.YLim(2)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
line([beta_poor beta_poor],[0 E_fig.Children.YLim(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig1_SD_LT_beta_range_softmax_overharvesting'],'-dsvg')


%% Panel: expected leave times for range of epsilon for E-GREEDY model

load([dataPath 'expectedLT_e-greedy_leheron'])

% beta prediction for MVT
ixRich = find(round(data.E_leave(:,2),1) >= round(optLT(1,2),1),1,"first");
ixPoor = find(round(data.E_leave(:,2),1) >= round(optLT(2,2),1),1,"first");
beta_rich = data.beta(ixRich);
beta_poor = data.beta(ixPoor);

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
colororder(color.patch);
semilogx(data.beta,data.E_leave,'LineWidth',widths.plot); hold on
xlabel('Epsilon (higher = explore)')
ylabel('Expected leaving time (s)')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig1_LT_beta_range_e_greedy_overharvesting'],'-dsvg')

%% Panel: expected e-greedy overharvesting leave times vs. optimal MVT 

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
plot(1:3,optLT(1,:),lines.exp,'Color',color.rich,'LineWidth',widths.plot); hold on
plot(1:3,optLT(2,:),lines.exp,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('Patch leaving time (s)')
set(gca,'XTick',1:3,'XTickLabel',{'Low', 'Medium', 'High'},'XLim',[0 4],'YLim',[0 30])

plot(1:3,data.E_leave(ixRich,:),lines.model,'Color',color.rich,'LineWidth',widths.plot); hold on
plot(1:3,data.E_leave(ixPoor,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig1_MVT_expectedLT_e-greedy_overharvesting'],'-dsvg')

%% Panel: expected leave times for range of epsilon for E-SOFTMAX (LAPSE) model

load([dataPath 'expectedLT_lapse_leheron_001'])

% beta prediction for MVT
ixRich = find(round(data.E_leave(:,2),1) >= round(optLT(1,2),1),1,"first");
ixPoor = find(round(data.E_leave(:,2),1) >= round(optLT(2,2),1),1,"first");
beta_rich = data.beta(ixRich);
beta_poor = data.beta(ixPoor);

E_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
colororder(color.patch);
semilogx(data.beta,data.E_leave,'LineWidth',widths.plot); hold on
xlabel('Beta (higher = exploit)')
ylabel('Expected leaving time (s)')

% plot overharvesting betas 
line([beta_rich beta_rich],[0 E_fig.Children.YLim(2)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
line([beta_poor beta_poor],[0 E_fig.Children.YLim(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig1_LT_beta_range_lapse_overharvesting_small_lapse'],'-dsvg')

load([dataPath 'expectedLT_lapse_leheron_01'])

% beta prediction for MVT
ixRich = find(round(data.E_leave(:,2),1) >= round(optLT(1,2),1),1,"first");
ixPoor = find(round(data.E_leave(:,2),1) >= round(optLT(2,2),1),1,"first");
beta_rich = data.beta(ixRich);
beta_poor = data.beta(ixPoor);

E_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
colororder(color.patch);
semilogx(data.beta,data.E_leave,'LineWidth',widths.plot); hold on
xlabel('Beta (higher = exploit)')
ylabel('Expected leaving time (s)')

% plot overharvesting betas 
line([beta_rich beta_rich],[0 E_fig.Children.YLim(2)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
line([beta_poor beta_poor],[0 E_fig.Children.YLim(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig1_LT_beta_range_lapse_overharvesting_mid_lapse'],'-dsvg')

load([dataPath 'expectedLT_lapse_leheron_04'])

% beta prediction for MVT
ixRich = find(round(data.E_leave(:,2),1) >= round(optLT(1,2),1),1,"first");
ixPoor = find(round(data.E_leave(:,2),1) >= round(optLT(2,2),1),1,"first");
beta_rich = data.beta(ixRich);
beta_poor = data.beta(ixPoor);

E_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
colororder(color.patch);
semilogx(data.beta,data.E_leave,'LineWidth',widths.plot); hold on
xlabel('Beta (higher = exploit)')
ylabel('Expected leaving time (s)')

% plot overharvesting betas 
line([beta_rich beta_rich],[0 E_fig.Children.YLim(2)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
line([beta_poor beta_poor],[0 E_fig.Children.YLim(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig1_LT_beta_range_lapse_overharvesting_large_lapse'],'-dsvg')


%% Panel: expected SD of leave times for range of epsilon for E-SOFTMAX (LAPSE) model

load([dataPath 'expectedLT_lapse_leheron_001'])

% beta prediction for MVT
ixRich = find(round(data.E_leave(:,2),1) >= round(optLT(1,2),1),1,"first");
ixPoor = find(round(data.E_leave(:,2),1) >= round(optLT(2,2),1),1,"first");
beta_rich = data.beta(ixRich);
beta_poor = data.beta(ixPoor);

SD_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
colororder(color.patch);
semilogx(data.beta,data.SD_leave,'LineWidth',widths.plot); hold on
xlabel('Beta (higher = exploit)')
ylabel('SD leaving times (s)')
ylim([0 15])

% plot overharvesting betas 
line([beta_rich beta_rich],[0 SD_fig.Children.YLim(2)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
line([beta_poor beta_poor],[0 SD_fig.Children.YLim(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig1_LT_beta_range_lapse_overharvesting_small_lapse_SD'],'-dsvg')

load([dataPath 'expectedLT_lapse_leheron_01'])

% beta prediction for MVT
ixRich = find(round(data.E_leave(:,2),1) >= round(optLT(1,2),1),1,"first");
ixPoor = find(round(data.E_leave(:,2),1) >= round(optLT(2,2),1),1,"first");
beta_rich = data.beta(ixRich);
beta_poor = data.beta(ixPoor);

SD_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
colororder(color.patch);
semilogx(data.beta,data.SD_leave,'LineWidth',widths.plot); hold on
xlabel('Beta (higher = exploit)')
ylabel('SD leaving times (s)')
ylim([0 15])

% plot overharvesting betas 
line([beta_rich beta_rich],[0 SD_fig.Children.YLim(2)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
line([beta_poor beta_poor],[0 SD_fig.Children.YLim(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig1_LT_beta_range_lapse_overharvesting_mid_lapse_SD'],'-dsvg')

load([dataPath 'expectedLT_lapse_leheron_04'])

% beta prediction for MVT
ixRich = find(round(data.E_leave(:,2),1) >= round(optLT(1,2),1),1,"first");
ixPoor = find(round(data.E_leave(:,2),1) >= round(optLT(2,2),1),1,"first");
beta_rich = data.beta(ixRich);
beta_poor = data.beta(ixPoor);

SD_fig = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
colororder(color.patch);
semilogx(data.beta,data.SD_leave,'LineWidth',widths.plot); hold on
xlabel('Beta (higher = exploit)')
ylabel('SD leaving times (s)')
ylim([0 15])

% plot overharvesting betas 
line([beta_rich beta_rich],[0 SD_fig.Children.YLim(2)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
line([beta_poor beta_poor],[0 SD_fig.Children.YLim(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig1_LT_beta_range_lapse_overharvesting_large_lapse_SD'],'-dsvg')


