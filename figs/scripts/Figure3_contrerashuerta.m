%%% Figure 3 foraging paper - Contreras-Huerta analysis
% Emma Scholey 7 Aug 2023 

clearvars; close all

dataPath = '../../data/fig_data/fig3/';

addpath('./plotFunctions')

run figure_properties_foraging.m

%% Panel: reward rates - exponential decay function

% reward rate values (from Le Heron) 
optimalAvgRR = [23.7388 19.2564]; 
r0=[34.5 57.5]; % initial yield 
a=0.11; % decay rate

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
print([export_path 'fig3_ContrerasHuerta_RR'],'-dsvg')
print([overleaf_path, 'fig3/', 'fig3_ContrerasHuerta_RR'],'-dsvg')

%% Panel: expected leave times for range of beta (softmax model) vs ContrerasHuerta data

load([dataPath 'expectedLT_softmax_contrerashuerta'])
load([dataPath 'subjectLT_contrerashuerta'])

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
print([export_path 'fig3_ContrerasHuerta_expectedLT'],'-dsvg')
print([overleaf_path, 'fig3/', 'fig3_ContrerasHuerta_expectedLT'],'-dsvg')

%% Panel: expected standard deviation for range of beta (softmax model) vs contrerashuerta data

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
print([export_path 'fig3_ContrerasHuerta_expectedSD'],'-dsvg')
print([overleaf_path, 'fig3/', 'fig3_ContrerasHuerta_expectedSD'],'-dsvg')


%% Panel: participant leave times vs fit softmax model

subjSEM = std(subject_leave_times.mean, [], 3)./sqrt(size(subject_leave_times.mean,3));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
errorbar(1:2,meanLT(1,:),subjSEM(1,:),lines.exp,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:2,meanLT(2,:),subjSEM(2,:),lines.exp,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('Patch leaving time (s)')
set(gca,'XTick',1:2,'XTickLabel',{'Low', 'High'},'XLim',[0 3],'YLim',[0 30])

% FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
% print([export_path 'fig3_ContrerasHuerta_subjectLT_softmax'],'-dsvg')

% add on MVT predictions to show overharvesting effect 
for e=1:2 %each env
    for p=1:2 % each patch type
        optLT(e,p)=(log(optimalAvgRR(e)/r0(p)))/-a; 
    end
end
plot(1:2,optLT(1,:),':','Color',color.rich,'LineWidth',widths.plot); 
plot(1:2,optLT(2,:),':','Color',color.poor,'LineWidth',widths.plot)


load([dataPath 'modelLT_contrerashuerta_separate_M1'])
meanLT_simulated = mean(subject_leave_times_simulated.mean,3);
simSEM = std(subject_leave_times_simulated.mean, [], 3)./sqrt(size(subject_leave_times_simulated.mean,3));

% figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
errorbar(1:2,meanLT_simulated(1,:),simSEM(1,:),lines.model,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:2,meanLT_simulated(2,:),simSEM(2,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('Patch leaving time (s)')
set(gca,'XTick',1:2,'XTickLabel',{'Low', 'High'},'XLim',[0 3],'YLim',[0 30])

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig3_ContrerasHuerta_simulatedLT_softmax'],'-dsvg')
print([overleaf_path, 'fig3/', 'fig3_ContrerasHuerta_simulatedLT_softmax'],'-dsvg')


%% Panel: BIC comparison for model
load([dataPath 'BIC_separate_contrerashuerta_choice_policy.mat'])
sumBIC = sum(models_BIC{:,1:2},2); % sum across environments

load([dataPath 'BIC_combined_contrerashuerta_fixedbeta.mat']) % append fixed beta model BIC
sumBIC = [sumBIC(1:2); models_BIC{:,1}]; 

modelNames = {'softmax', 'greedy', 'softmax fixed'};

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
bar(sumBIC,'FaceColor',color.general, 'EdgeColor', color.general); hold on

% find best model and highlight
minBIC = min(sumBIC);
bestModel = find(sumBIC == minBIC);
bar(bestModel,minBIC, 'FaceColor',color.highlight, 'EdgeColor', color.highlight)

ylabel('BIC (sum)')
set(gca,'XTickLabel',modelNames, 'YLim', [18000 25000], 'XTickLabelRotation',50)

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig3_ContrerasHuerta_BIC'],'-dsvg')
print([overleaf_path, 'fig3/', 'fig3_ContrerasHuerta_BIC'],'-dsvg')


%% Panel: distribution of beta fit with LT for each environment
load([dataPath 'fitting_results_separate_M1_contrerashuerta.mat'])

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square); hold on

plot individual points
scatter(ones(size(minNLLFitParams_rich.beta)), minNLLFitParams_rich.beta,'o','MarkerFaceColor', color.rich,'MarkerEdgeColor', 'w'); 
scatter(2*ones(size(minNLLFitParams_poor.beta)), minNLLFitParams_poor.beta,'o','MarkerFaceColor', color.poor,'MarkerEdgeColor', 'w'); 

plot error bars 
betaSEM_rich = std(minNLLFitParams_rich.beta)./sqrt(numel(minNLLFitParams_rich.beta));
betaSEM_poor = std(minNLLFitParams_poor.beta)./sqrt(numel(minNLLFitParams_poor.beta));

errorbar(1,mean(minNLLFitParams_rich.beta),betaSEM_rich, '_', 'LineWidth',widths.plot,'Color',[0.1 0.1 0.1],'CapSize', 0); 
errorbar(2,mean(minNLLFitParams_poor.beta),betaSEM_poor, '_','LineWidth',widths.plot,'Color',[0.1 0.1 0.1],'CapSize', 0); 

set(gca,'XTick',1:2, 'XTickLabel', {'Rich', 'Poor'}, 'XLim', [0,3] ,'YScale', 'log', 'YLim', [0,1])
ylabel('Fit beta (higher = exploit)')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig3_ContrerasHuerta_fit_beta'],'-dsvg')
print([overleaf_path, 'fig3/', 'fig3_ContrerasHuerta_fit_beta'],'-dsvg')

%% Panel: correlation of beta fit with LT for each environment

% meanLT_rich = mean(squeeze(subject_leave_times.mean(1,:,:)));
% meanLT_poor = mean(squeeze(subject_leave_times.mean(2,:,:)));
% 
% figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
% scatter(minNLLFitParams_rich.beta, meanLT_rich,'MarkerFaceColor',color.rich,'MarkerEdgeColor', color.rich, 'MarkerFaceAlpha',M_alpha, 'MarkerEdgeAlpha',M_alpha); hold on
% scatter(minNLLFitParams_poor.beta, meanLT_poor,'MarkerFaceColor',color.poor,'MarkerEdgeColor', color.poor, 'MarkerFaceAlpha',M_alpha, 'MarkerEdgeAlpha',M_alpha); hold on
% set(gca,'XLim', [0 2],'XScale', 'log','XTick', [10^-2, 10^-1, 10^0, 10^1])
% xlabel('Fit beta (higher = exploit)')
% ylabel('Patch leaving time (s)')
% 
% FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
% print([export_path 'fig3_ContrerasHuerta_fit_beta_LT_correlation'],'-dsvg')
% 
% % statistics
% [pearson_coeff, p_val]  = corr(log([minNLLFitParams_rich.beta; minNLLFitParams_poor.beta]), [meanLT_rich'; meanLT_poor'])
% 
% [H0_reject, p_val_t_test, ~, stats] = ttest(log(minNLLFitParams_rich.beta), log(minNLLFitParams_poor.beta))
% print([export_path 'fig3_ContrerasHuerta_fit_beta_LT_correlation'],'-dsvg')

%% Panel: Participant SD vs fit softmax model 

meanLT_SD_simulated = mean(subject_leave_times_simulated.sd,3);

subjSE_SD = std(subject_leave_times.sd, [], 3)./sqrt(size(subject_leave_times.sd,3));
simSE_SD = std(subject_leave_times_simulated.sd, [], 3)./sqrt(size(subject_leave_times_simulated.sd,3));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
errorbar(1:2,meanSD(1,:),subjSE_SD(1,:),lines.exp,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:2,meanSD(2,:),subjSE_SD(2,:),lines.exp,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('SD leaving time (s)')
set(gca,'XTick',1:2,'XTickLabel',{'Low', 'High'},'XLim',[0 3],'YLim',[0 10])

% % FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
% print([export_path 'fig3_ContrerasHuerta_subjectLT_SD'],'-dsvg')

% figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
errorbar(1:2,meanLT_SD_simulated(1,:),simSE_SD(1,:),lines.model,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:2,meanLT_SD_simulated(2,:),simSE_SD(2,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('SD leaving time (s)')
set(gca,'XTick',1:2,'XTickLabel',{'Low', 'High'},'XLim',[0 3],'YLim',[0 10])

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig3_ContrerasHuerta_simulatedLT_SD'],'-dsvg')
print([overleaf_path, 'fig3/', 'fig3_ContrerasHuerta_simulatedLT_SD_softmax'],'-dsvg')


%% SUPPLEMENTARY FIGURES %%%%%%%%%%%%%%%%%%%%%%%%
%% Panel: participant leave times vs fit softmax model

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
% participant data
errorbar(1:2,meanLT(1,:),subjSEM(1,:),lines.exp,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:2,meanLT(2,:),subjSEM(2,:),lines.exp,'Color',color.poor,'LineWidth',widths.plot)

load([dataPath 'modelLT_contrerashuerta_combined_M1'])
% bad model data
meanLT_simulated = mean(subject_leave_times_simulated.mean,3);
simSEM = std(subject_leave_times_simulated.mean, [], 3)./sqrt(size(subject_leave_times_simulated.mean,3));

errorbar(1:2,meanLT_simulated(1,:),simSEM(1,:),lines.model,'Color',color.rich,'LineWidth',widths.plot)
errorbar(1:2,meanLT_simulated(2,:),simSEM(2,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)

xlabel('Patch yield')
ylabel('Patch leaving time (s)')
set(gca,'XTick',1:2,'XTickLabel',{'Low', 'High'},'XLim',[0 3],'YLim',[0 30])

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig3_ContrerasHuerta_simulatedLT_softmax_fixedbeta'],'-dsvg')
print([overleaf_path, 'fig3/', 'fig3_ContrerasHuerta_simulatedLT_softmax_fixedbeta'],'-dsvg')

%% Panel: lapse fits are near to 0 and do not correlate with leaving times 
% load([dataPath 'fitting_results_separate_M4_contrerashuerta.mat'])
% 
% figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
% scatter(minNLLFitParams_rich.epsilon, meanLT_rich,'MarkerFaceColor',color.rich,'MarkerEdgeColor', color.rich, 'MarkerFaceAlpha',M_alpha, 'MarkerEdgeAlpha',M_alpha); hold on
% scatter(minNLLFitParams_poor.epsilon, meanLT_poor,'MarkerFaceColor',color.poor,'MarkerEdgeColor', color.poor, 'MarkerFaceAlpha',M_alpha, 'MarkerEdgeAlpha',M_alpha); hold on
% set(gca,'XScale', 'log')
% xlabel('Fit epsilon (higher = more lapses)')
% ylabel('Patch leaving time (s)')
% 
% FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
% print([export_path 'fig3_ContrerasHuerta_fit_lapse_LT_correlation'],'-dsvg')
%print([overleaf_path, 'fig3/', 'fig3_ContrerasHuerta_fit_lapse_LT_correlation'],'-dsvg')

%% STATISTICS %%%%%%%%%%%%%%%%%
% statistics - mean LT - confirms dataset results
rich_mean = squeeze(subject_leave_times.mean(1,:,:))';
poor_mean = squeeze(subject_leave_times.mean(2,:,:))';

mean_table = [rich_mean;poor_mean];

[p_val, tbl] = anova2(mean_table,size(rich_mean,1)); % no significant difference between patch/environment for LT 

% statistics - SD LT - confirms novel SD predictions 
rich_sd = squeeze(subject_leave_times.sd(1,:,:))';
poor_sd = squeeze(subject_leave_times.sd(2,:,:))';

sd_table = [rich_sd;poor_sd];

[p_val, tbl] = anova2(sd_table,size(rich_sd,1)); % no significant difference between patch/environment for LT 

% statistics - SD LT SIMULATIONS
rich_sd = squeeze(subject_leave_times_simulated.sd(1,:,:))';
poor_sd = squeeze(subject_leave_times_simulated.sd(2,:,:))';

sd_table = [rich_sd;poor_sd];

[p_val, tbl] = anova2(sd_table,size(rich_sd,1)); 
