%%% Figure 2 foraging paper - Le Heron analysis
% Emma Scholey 7 Aug 2023 

clearvars; close all

dataPath = '../../data/fig_data/fig2/';

addpath('./plotFunctions')

run figure_properties_foraging.m

%% Panel: reward rates - exponential decay function

% reward rate values (from Le Heron) 
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

line([0 RR_fig.Children.XLim(2)],[optimalAvgRR(1) optimalAvgRR(1)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
line([0 RR_fig.Children.XLim(2)],[optimalAvgRR(2) optimalAvgRR(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig2_LeHeron_RR'],'-dsvg')

%% Panel: expected leave times for range of beta (softmax model) vs leheron data

load([dataPath 'expectedLT_softmax_leheron'])
load([dataPath 'subjectLT_leheron'])

meanLT = mean(subject_leave_times.mean,3);
meanSD = mean(subject_leave_times.sd,3);

% beta prediction for participant group mean
ixRich = find(round(data.E_leave(:,2),1) >= round(meanLT(1,2),1),1,"first");
ixPoor = find(round(data.E_leave(:,2),1) >= round(meanLT(2,2),1),1,"first");

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
print([export_path 'fig2_LeHeron_expectedLT'],'-dsvg')

%% Panel: expected standard deviation for range of beta (softmax model) vs leheron data

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
colororder(color.patch);
semilogx(data.beta,data.SD_leave,'LineWidth',widths.plot); hold on
xlabel('Beta (higher = exploit)')
ylabel('SD leaving times (s)')
ylim([0 15])
% plot overharvesting betas 
line([beta_rich beta_rich],[0 E_fig.Children.YLim(2)],'Color',color.rich, 'LineWidth',widths.plot, 'LineStyle', '--')
line([beta_poor beta_poor],[0 E_fig.Children.YLim(2)],'Color',color.poor, 'LineWidth',widths.plot, 'LineStyle', '--')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig2_LeHeron_expectedSD'],'-dsvg')

%% Panel: participant leave times vs fit softmax model

subjSEM = std(subject_leave_times.mean, [], 3)./sqrt(size(subject_leave_times.mean,3));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
errorbar(1:3,meanLT(1,:),subjSEM(1,:),lines.exp,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:3,meanLT(2,:),subjSEM(2,:),lines.exp,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('Patch leaving time (s)')
set(gca,'XTick',1:3,'XTickLabel',{'Low', 'Medium', 'High'},'XLim',[0 4],'YLim',[0 30])

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig2_LeHeron_subjectLT_softmax'],'-dsvg')


load([dataPath 'modelLT_leheron_M1'])
meanLT_simulated = mean(subject_leave_times_simulated.mean,3);
simSEM = std(subject_leave_times_simulated.mean, [], 3)./sqrt(size(subject_leave_times_simulated.mean,3));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
errorbar(1:3,meanLT_simulated(1,:),simSEM(1,:),lines.model,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:3,meanLT_simulated(2,:),simSEM(2,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('Patch leaving time (s)')
set(gca,'XTick',1:3,'XTickLabel',{'Low', 'Medium', 'High'},'XLim',[0 4],'YLim',[0 30])

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig2_LeHeron_simulatedLT_softmax'],'-dsvg')


%% Panel: BIC comparison for model
load([dataPath 'BIC_separate_leheron.mat'])

meanBIC = mean(models_BIC{:,1:2},2); % take average across environments if separate fitting
modelNames = {'softmax', 'e-greedy','lapse', 'RW AvgRR', 'RW Patch+AvgRR', 'RW PatchRR'};

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
bar(meanBIC,'FaceColor',color.general, 'EdgeColor', color.general); hold on

% find best model and highlight
minBIC = min(meanBIC);
bestModel = find(meanBIC == minBIC);
bar(bestModel,minBIC, 'FaceColor',color.highlight, 'EdgeColor', color.highlight)

ylabel('BIC (sum)')
set(gca,'XTickLabel',modelNames, 'YLim', [5500 8000], 'XTickLabelRotation',50)

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig2_LeHeron_BIC'],'-dsvg')

%% Panel: distribution of beta fit with LT for each environment
load([dataPath 'fitting_results_separate_M1_leheron.mat'])

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
bar(1,mean(minNLLFitParams_rich.beta),'FaceColor',color.rich, 'EdgeColor',color.rich); hold on 
bar(2,mean(minNLLFitParams_poor.beta),'FaceColor',color.poor, 'EdgeColor',color.poor); 

% plot individual points
scatter(ones(size(minNLLFitParams_rich.beta)), minNLLFitParams_rich.beta,'o','MarkerFaceColor', marker.rich,'MarkerEdgeColor', marker.rich,'MarkerEdgeAlpha', M_alpha,'MarkerFaceAlpha', M_alpha); 
scatter(2*ones(size(minNLLFitParams_poor.beta)), minNLLFitParams_poor.beta,'o','MarkerFaceColor', marker.poor,'MarkerEdgeColor', marker.poor,'MarkerEdgeAlpha', M_alpha,'MarkerFaceAlpha', M_alpha); 

% plot error bars 
betaSEM_rich = std(minNLLFitParams_rich.beta)./sqrt(numel(minNLLFitParams_rich.beta));
betaSEM_poor = std(minNLLFitParams_poor.beta)./sqrt(numel(minNLLFitParams_poor.beta));

errorbar(1,mean(minNLLFitParams_rich.beta),betaSEM_rich,'LineStyle', 'none','LineWidth',widths.plot,'Color',[0.1 0.1 0.1],'CapSize', 0); 
errorbar(2,mean(minNLLFitParams_poor.beta),betaSEM_poor,'LineStyle', 'none','LineWidth',widths.plot,'Color',[0.1 0.1 0.1],'CapSize', 0); 


[H0_reject, p_val_t_test] = ttest(log(minNLLFitParams_rich.beta), log(minNLLFitParams_poor.beta))

set(gca,'XTick',1:2, 'XTickLabel', {'Rich', 'Poor'}, 'YScale', 'log')
ylabel('Fit beta (higher = exploit)')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig2_LeHeron_fit_beta_LT_group'],'-dsvg')



%% Panel: correlation of beta fit with LT for each environment

meanLT_rich = mean(squeeze(subject_leave_times.mean(1,:,:)));
meanLT_poor = mean(squeeze(subject_leave_times.mean(2,:,:)));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
scatter(minNLLFitParams_rich.beta, meanLT_rich,'MarkerFaceColor',color.rich,'MarkerEdgeColor', color.rich, 'MarkerFaceAlpha',M_alpha, 'MarkerEdgeAlpha',M_alpha); hold on
scatter(minNLLFitParams_poor.beta, meanLT_poor,'MarkerFaceColor',color.poor,'MarkerEdgeColor', color.poor, 'MarkerFaceAlpha',M_alpha, 'MarkerEdgeAlpha',M_alpha); hold on
set(gca,'XLim', [0 2],'XScale', 'log','XTick', [10^-2, 10^-1, 10^0, 10^1])
xlabel('Fit beta (higher = exploit)')
ylabel('Patch leaving time (s)')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig2_LeHeron_fit_beta_LT_correlation'],'-dsvg')

% statistics
[pearson_coeff_rich, p_val_rich]  = corr(minNLLFitParams_rich.beta, meanLT_rich')
[pearson_coeff_poor, p_val_poor]  = corr(minNLLFitParams_poor.beta, meanLT_poor')

%% Panel: Participant SD vs fit softmax model 

meanLT_SD_simulated = mean(subject_leave_times_simulated.sd,3);

subjSE_SD = std(subject_leave_times.sd, [], 3)./sqrt(size(subject_leave_times.sd,3));
simSE_SD = std(subject_leave_times_simulated.sd, [], 3)./sqrt(size(subject_leave_times_simulated.sd,3));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
errorbar(1:3,meanSD(1,:),subjSE_SD(1,:),lines.exp,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:3,meanSD(2,:),subjSE_SD(2,:),lines.exp,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('SD leaving times (s)')
set(gca,'XTick',1:3,'XTickLabel',{'Low', 'Medium', 'High'},'XLim',[0 4],'YLim',[0 10])

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig2_LeHeron_subjectLT_SD_softmax'],'-dsvg')

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
errorbar(1:3,meanLT_SD_simulated(1,:),simSE_SD(1,:),lines.model,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:3,meanLT_SD_simulated(2,:),simSE_SD(2,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('SD leaving times (s)')
set(gca,'XTick',1:3,'XTickLabel',{'Low', 'Medium', 'High'},'XLim',[0 4],'YLim',[0 10])

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig2_LeHeron_simulatedLT_SD_softmax'],'-dsvg')

% statistics - mean LT - confirms le heron results
rich_mean = squeeze(subject_leave_times.mean(1,:,:))';
poor_mean = squeeze(subject_leave_times.mean(2,:,:))';

mean_table = [rich_mean;poor_mean];

[p_val, tbl] = anova2(mean_table,size(rich_mean,1)); 

% statistics - SD LT - confirms novel SD predictions 
rich_sd = squeeze(subject_leave_times.sd(1,:,:))';
poor_sd = squeeze(subject_leave_times.sd(2,:,:))';

sd_table = [rich_sd;poor_sd];

[p_val, tbl] = anova2(sd_table,size(rich_sd,1)); 

% statistics - SD LT SIMULATIONS
rich_sd = squeeze(subject_leave_times_simulated.sd(1,:,:))';
poor_sd = squeeze(subject_leave_times_simulated.sd(2,:,:))';

sd_table = [rich_sd;poor_sd];

[p_val, tbl] = anova2(sd_table,size(rich_sd,1)); 