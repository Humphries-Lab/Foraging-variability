%%% Figure 4 foraging paper - beta function analysis
% Emma Scholey 10 Aug 2023 

clearvars; close all

dataPath = '../../data/fig_data/fig5/';

addpath('./plotFunctions')

run figure_properties_foraging.m

%% Panel: BIC comparison for model
load([dataPath 'BIC_combined_leheron_functions.mat'])

modelNames = {'scalar MVT', 'scalar RR',...
              'exp MVT', 'exp RR',...
              'hyperbolic MVT', 'hyperbolic RR',...
              '2 scalar MVT', '2 scalar RR',...
              '2 exp MVT', '2 exp RR'};

h = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.rectangle);
bar(models_BIC.combined(1:10,:),'FaceColor',color.general, 'EdgeColor', color.general); hold on

% find best model and highlight
minBIC = min(models_BIC.combined);
bestModel = find(models_BIC.combined == minBIC);
bar(bestModel,minBIC, 'FaceColor',color.highlight, 'EdgeColor', color.highlight)
ylabel('BIC (sum)')
set(gca,'XTick',[1:10],'XTickLabel',modelNames, 'YLim', [12000, 13000],'XTickLabelRotation',90)

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_beta_function_BIC'],'-dsvg')

%% Panel: participant leave times vs best fit beta function model

load([dataPath 'subjectLT_leheron'])
meanLT = mean(subject_leave_times.mean,3);
meanSD = mean(subject_leave_times.sd,3);
subjSEM = std(subject_leave_times.mean, [], 3)./sqrt(size(subject_leave_times.mean,3));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
errorbar(1:3,meanLT(1,:),subjSEM(1,:),lines.exp,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:3,meanLT(2,:),subjSEM(2,:),lines.exp,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('Patch leaving time (s)')
set(gca,'XTick',1:3,'XTickLabel',{'Low', 'Medium', 'High'},'XLim',[0 4],'YLim',[0 30])

% put on same plot 
load([dataPath 'modelLT_leheron_combined_M9'])
meanLT_simulated = mean(subject_leave_times_simulated.mean,3);
simSEM = std(subject_leave_times_simulated.mean, [], 3)./sqrt(size(subject_leave_times_simulated.mean,3));

errorbar(1:3,meanLT_simulated(1,:),simSEM(1,:),lines.model,'Color',color.rich,'LineWidth',widths.plot);
errorbar(1:3,meanLT_simulated(2,:),simSEM(2,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_LeHeron_simulatedLT_beta_function'],'-dsvg')

%% Panel: correlation of lambda fit with LT for each environment
load([dataPath 'fitting_results_combined_M9_leheron.mat'])

meanLT_subj = squeeze(mean(mean(subject_leave_times.mean(:,:,:),1),2));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
scatter(minNLLFitParams.lambda, meanLT_subj,'MarkerFaceColor',color.general,'MarkerEdgeColor', color.general, 'MarkerFaceAlpha',M_alpha, 'MarkerEdgeAlpha',M_alpha); hold on
set(gca,'XScale', 'log','YLim', [0 50])
xlabel('Fit lambda (higher = exploit)')
ylabel('Patch leaving time (s)')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_LeHeron_fit_lambda_LT_correlation'],'-dsvg')

% statistics
[pearson_coeff, p_val]  = corr(log(minNLLFitParams.lambda), meanLT_subj)

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

errorbar(1:3,meanLT_SD_simulated(1,:),simSE_SD(1,:),lines.model,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:3,meanLT_SD_simulated(2,:),simSE_SD(2,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_LeHeron_simulatedLT_SD_beta_function'],'-dsvg')

% % statistics - SD LT - confirms novel SD predictions 
% rich_sd = squeeze(subject_leave_times.sd(1,:,:))';
% poor_sd = squeeze(subject_leave_times.sd(2,:,:))';
% 
% sd_table = [rich_sd;poor_sd];
% 
% [p_val, tbl] = anova2(sd_table,size(rich_sd,1)); 
% 
% % statistics - SD LT SIMULATIONS
% rich_sd = squeeze(subject_leave_times_simulated.sd(1,:,:))';
% poor_sd = squeeze(subject_leave_times_simulated.sd(2,:,:))';
% 
% sd_table = [rich_sd;poor_sd];
% 
% [p_val, tbl] = anova2(sd_table,size(rich_sd,1)); 


%% SUPPLEMENTARY Panel: BIC comparison for contreras-huerta datasest 
load([dataPath 'BIC_combined_contrerashuerta_functions.mat'])

h = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.rectangle);
bar(models_BIC.combined(1:10,:),'FaceColor',color.general, 'EdgeColor', color.general); hold on

% find best model and highlight
minBIC = min(models_BIC.combined);
bestModel = find(models_BIC.combined == minBIC);
bar(bestModel,minBIC, 'FaceColor',color.highlight, 'EdgeColor', color.highlight)
ylabel('BIC (sum)')
set(gca,'XTick',[1:10],'XTickLabel',modelNames, 'YLim',[18500, 19500],'XTickLabelRotation',90)

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_beta_function_BIC_contrerashuerta'],'-dsvg')

%% SUPPLEMENTARY Panel: BIC comparison for kane datasest 
load([dataPath 'BIC_combined_kane_functions.mat'])

h = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.rectangle);
bar(models_BIC.combined(1:10,:),'FaceColor',color.general, 'EdgeColor', color.general); hold on

% find best model and highlight
minBIC = min(models_BIC.combined);
bestModel = find(models_BIC.combined == minBIC);
bar(bestModel,minBIC, 'FaceColor',color.highlight, 'EdgeColor', color.highlight)
ylabel('BIC (sum)')
set(gca,'XTick',[1:10],'XTickLabel',modelNames, 'YLim',[11000, 11500],'XTickLabelRotation',90)

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_beta_function_BIC_kane'],'-dsvg')

%% SUPPLEMENTARY Panel: contrerashuerta participant leave times vs best fit beta function model
load([dataPath 'modelLT_contrerashuerta_combined_M9'])
meanLT_simulated = mean(subject_leave_times_simulated.mean,3);
simSEM = std(subject_leave_times_simulated.mean, [], 3)./sqrt(size(subject_leave_times_simulated.mean,3));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
errorbar(1:2,meanLT_simulated(1,:),simSEM(1,:),lines.model,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:2,meanLT_simulated(2,:),simSEM(2,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('Patch leaving time (s)')
set(gca,'XTick',1:2,'XTickLabel',{'Low', 'High'},'XLim',[0 3],'YLim',[0 30])

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_ContrerasHuerta_simulatedLT_beta_function'],'-dsvg')


% Look at variance too 
meanLT_SD_simulated = mean(subject_leave_times_simulated.sd,3);
simSE_SD = std(subject_leave_times_simulated.sd, [], 3)./sqrt(size(subject_leave_times_simulated.sd,3));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
errorbar(1:2,meanLT_SD_simulated(1,:),simSE_SD(1,:),lines.model,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:2,meanLT_SD_simulated(2,:),simSE_SD(2,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('SD leaving times (s)')
set(gca,'XTick',1:2,'XTickLabel',{'Low', 'High'},'XLim',[0 3],'YLim',[0 10])

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_ContrerasHuerta_simulatedLT_SD_beta_function'],'-dsvg')


%% SUPPLEMENTARY Panel: kane leave times vs best fit beta function model
load([dataPath 'modelLT_kane_combined_M9'])
meanLT_simulated = mean(subject_leave_times_simulated.mean,3);
simSEM = std(subject_leave_times_simulated.mean, [], 3)./sqrt(size(subject_leave_times_simulated.mean,3));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
errorbar(1:3,meanLT_simulated(1,:),simSEM(1,:),lines.model,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:3,meanLT_simulated(2,:),simSEM(2,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('Patch leaving time (s)')
set(gca,'XTick',1:3,'XTickLabel',{'Low', 'Medium', 'High'},'XLim',[0 4],'YLim',[0 30])

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_Kane_simulatedLT_beta_function'],'-dsvg')

% Look at variance too 
meanLT_SD_simulated = mean(subject_leave_times_simulated.sd,3);
simSE_SD = std(subject_leave_times_simulated.sd, [], 3)./sqrt(size(subject_leave_times_simulated.sd,3));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
errorbar(1:3,meanLT_SD_simulated(1,:),simSE_SD(1,:),lines.model,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:3,meanLT_SD_simulated(2,:),simSE_SD(2,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('SD leaving times (s)')
set(gca,'XTick',1:3,'XTickLabel',{'Low', 'Medium', 'High'},'XLim',[0 4],'YLim',[0 5])

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_Kane_simulatedLT_SD_beta_function'],'-dsvg')

%% SUPPLEMENTARY Panel: BIC comparison for Rescorla Wagner models 

%%%%% Le Heron 
load([dataPath 'BIC_combined_leheron_learning.mat'])
modelNames = {'scalar RW', 'exp RW', 'hyperbolic RW', '2 scalar RW', '2 exp RW'};

h = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
bar(models_BIC.combined(1:5,:),'FaceColor',color.general, 'EdgeColor', color.general); hold on

% find best model and highlight
minBIC = min(models_BIC.combined);
bestModel = find(models_BIC.combined == minBIC);
bar(bestModel,minBIC, 'FaceColor',color.highlight, 'EdgeColor', color.highlight)
ylabel('BIC (sum)')
set(gca,'XTick',[1:5],'XTickLabel',modelNames, 'YLim', [12000, 13000],'XTickLabelRotation',90)

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_leheron_beta_function_BIC_learning'],'-dsvg')

%%%%% Contreras Huerta
load([dataPath 'BIC_combined_contrerashuerta_learning.mat'])

h = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
bar(models_BIC.combined(1:5,:),'FaceColor',color.general, 'EdgeColor', color.general); hold on

% find best model and highlight
minBIC = min(models_BIC.combined);
bestModel = find(models_BIC.combined == minBIC);
bar(bestModel,minBIC, 'FaceColor',color.highlight, 'EdgeColor', color.highlight)
ylabel('BIC (sum)')
set(gca,'XTick',[1:5],'XTickLabel',modelNames, 'YLim', [18500, 19500],'XTickLabelRotation',90)

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_contrerashuerta_beta_function_BIC_learning'],'-dsvg')

%%%%% Contreras Huerta
load([dataPath 'BIC_combined_kane_learning.mat'])

h = figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
bar(models_BIC.combined(1:5,:),'FaceColor',color.general, 'EdgeColor', color.general); hold on

% find best model and highlight
minBIC = min(models_BIC.combined);
bestModel = find(models_BIC.combined == minBIC);
bar(bestModel,minBIC, 'FaceColor',color.highlight, 'EdgeColor', color.highlight)
ylabel('BIC (sum)')
set(gca,'XTick',[1:5],'XTickLabel',modelNames, 'YLim', [11250, 11750],'XTickLabelRotation',90)

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_kane_beta_function_BIC_learning'],'-dsvg')

%% SUPPLEMENTARY Panel: leave times for Rescorla Wagner models 

%%%%% Le Heron 
load([dataPath 'subjectLT_leheron'])
meanLT = mean(subject_leave_times.mean,3);
meanSD = mean(subject_leave_times.sd,3);
subjSEM = std(subject_leave_times.mean, [], 3)./sqrt(size(subject_leave_times.mean,3));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
errorbar(1:3,meanLT(1,:),subjSEM(1,:),lines.exp,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:3,meanLT(2,:),subjSEM(2,:),lines.exp,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('Patch leaving time (s)')
set(gca,'XTick',1:3,'XTickLabel',{'Low', 'Medium', 'High'},'XLim',[0 4],'YLim',[0 30])

% put on same plot 
load([dataPath 'modelLT_leheron_combined_M10'])
meanLT_simulated = mean(subject_leave_times_simulated.mean,3);
simSEM = std(subject_leave_times_simulated.mean, [], 3)./sqrt(size(subject_leave_times_simulated.mean,3));

errorbar(1:3,meanLT_simulated(1,:),simSEM(1,:),lines.model,'Color',color.rich,'LineWidth',widths.plot);
errorbar(1:3,meanLT_simulated(2,:),simSEM(2,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_LeHeron_simulatedLT_learning_function'],'-dsvg')

% do same for variance
meanLT_SD_simulated = mean(subject_leave_times_simulated.sd,3);

subjSE_SD = std(subject_leave_times.sd, [], 3)./sqrt(size(subject_leave_times.sd,3));
simSE_SD = std(subject_leave_times_simulated.sd, [], 3)./sqrt(size(subject_leave_times_simulated.sd,3));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
errorbar(1:3,meanSD(1,:),subjSE_SD(1,:),lines.exp,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:3,meanSD(2,:),subjSE_SD(2,:),lines.exp,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('SD leaving times (s)')
set(gca,'XTick',1:3,'XTickLabel',{'Low', 'Medium', 'High'},'XLim',[0 4],'YLim',[0 10])

errorbar(1:3,meanLT_SD_simulated(1,:),simSE_SD(1,:),lines.model,'Color',color.rich,'LineWidth',widths.plot)
errorbar(1:3,meanLT_SD_simulated(2,:),simSE_SD(2,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_LeHeron_simulatedLT_SD_learning_function'],'-dsvg')

%%%%% Contreras Huerta
load([dataPath 'subjectLT_contrerashuerta'])
meanLT = mean(subject_leave_times.mean,3);
meanSD = mean(subject_leave_times.sd,3);
subjSEM = std(subject_leave_times.mean, [], 3)./sqrt(size(subject_leave_times.mean,3));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
errorbar(1:2,meanLT(1,:),subjSEM(1,:),lines.exp,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:2,meanLT(2,:),subjSEM(2,:),lines.exp,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('Patch leaving time (s)')
set(gca,'XTick',1:2,'XTickLabel',{'Low', 'High'},'XLim',[0 3],'YLim',[0 30])

% put on same plot 
load([dataPath 'modelLT_contrerashuerta_combined_M10'])
meanLT_simulated = mean(subject_leave_times_simulated.mean,3);
simSEM = std(subject_leave_times_simulated.mean, [], 3)./sqrt(size(subject_leave_times_simulated.mean,3));

errorbar(1:2,meanLT_simulated(1,:),simSEM(1,:),lines.model,'Color',color.rich,'LineWidth',widths.plot);
errorbar(1:2,meanLT_simulated(2,:),simSEM(2,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_ContrerasHuerta_simulatedLT_learning_function'],'-dsvg')

% do same for variance
meanLT_SD_simulated = mean(subject_leave_times_simulated.sd,3);

subjSE_SD = std(subject_leave_times.sd, [], 3)./sqrt(size(subject_leave_times.sd,3));
simSE_SD = std(subject_leave_times_simulated.sd, [], 3)./sqrt(size(subject_leave_times_simulated.sd,3));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
errorbar(1:2,meanSD(1,:),subjSE_SD(1,:),lines.exp,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:2,meanSD(2,:),subjSE_SD(2,:),lines.exp,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('SD leaving times (s)')
set(gca,'XTick',1:2,'XTickLabel',{'Low', 'High'},'XLim',[0 3],'YLim',[0 10])

errorbar(1:2,meanLT_SD_simulated(1,:),simSE_SD(1,:),lines.model,'Color',color.rich,'LineWidth',widths.plot)
errorbar(1:2,meanLT_SD_simulated(2,:),simSE_SD(2,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_ContrerasHuerta_simulatedLT_SD_learning_function'],'-dsvg')

%%%%% Kane 
load([dataPath 'subjectLT_kane'])
meanLT = mean(subject_leave_times.mean,3);
meanSD = mean(subject_leave_times.sd,3);
subjSEM = std(subject_leave_times.mean, [], 3)./sqrt(size(subject_leave_times.mean,3));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
errorbar(1:3,meanLT(1,:),subjSEM(1,:),lines.exp,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:3,meanLT(2,:),subjSEM(2,:),lines.exp,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('Patch leaving time (s)')
set(gca,'XTick',1:3,'XTickLabel',{'Low', 'Medium', 'High'},'XLim',[0 4],'YLim',[0 30])

% put on same plot 
load([dataPath 'modelLT_kane_combined_M10'])
meanLT_simulated = mean(subject_leave_times_simulated.mean,3);
simSEM = std(subject_leave_times_simulated.mean, [], 3)./sqrt(size(subject_leave_times_simulated.mean,3));

errorbar(1:3,meanLT_simulated(1,:),simSEM(1,:),lines.model,'Color',color.rich,'LineWidth',widths.plot);
errorbar(1:3,meanLT_simulated(2,:),simSEM(2,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_Kane_simulatedLT_learning_function'],'-dsvg')

% do same for variance
meanLT_SD_simulated = mean(subject_leave_times_simulated.sd,3);

subjSE_SD = std(subject_leave_times.sd, [], 3)./sqrt(size(subject_leave_times.sd,3));
simSE_SD = std(subject_leave_times_simulated.sd, [], 3)./sqrt(size(subject_leave_times_simulated.sd,3));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
errorbar(1:3,meanSD(1,:),subjSE_SD(1,:),lines.exp,'Color',color.rich,'LineWidth',widths.plot); hold on
errorbar(1:3,meanSD(2,:),subjSE_SD(2,:),lines.exp,'Color',color.poor,'LineWidth',widths.plot)
xlabel('Patch yield')
ylabel('SD leaving times (s)')
set(gca,'XTick',1:3,'XTickLabel',{'Low', 'Medium', 'High'},'XLim',[0 4],'YLim',[0 5])

errorbar(1:3,meanLT_SD_simulated(1,:),simSE_SD(1,:),lines.model,'Color',color.rich,'LineWidth',widths.plot)
errorbar(1:3,meanLT_SD_simulated(2,:),simSE_SD(2,:),lines.model,'Color',color.poor,'LineWidth',widths.plot)

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_Kane_simulatedLT_SD_learning_function'],'-dsvg')

%% SUPPLEMENTARY Panel: correlation of learning rate fits vs leaving times for Rescorla Wagner models 

%%%%% Le Heron 
load([dataPath 'subjectLT_leheron'])
load([dataPath 'fitting_results_combined_M10_leheron.mat'])

meanLT_subj = squeeze(mean(mean(subject_leave_times.mean(:,:,:),1),2));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
scatter(minNLLFitParams.alphaRho, meanLT_subj,'MarkerFaceColor',color.general,'MarkerEdgeColor', color.general, 'MarkerFaceAlpha',M_alpha, 'MarkerEdgeAlpha',M_alpha); hold on
set(gca,'YLim', [0 50]) % small x axis, say in range of previous studies e.g. Garett & Daw. Check Lloyds studies too 
xlabel('Fit alpha')
ylabel('Patch leaving time (s)')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_LeHeron_fit_alpha_LT_correlation'],'-dsvg')

% statistics
[pearson_coeff, p_val]  = corr(minNLLFitParams.alphaRho, meanLT_subj)

%%%%% Contreras Huerta 
load([dataPath 'subjectLT_contrerashuerta'])
load([dataPath 'fitting_results_combined_M10_contrerashuerta.mat'])

meanLT_subj = squeeze(mean(mean(subject_leave_times.mean(:,:,:),1),2));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
scatter(minNLLFitParams.alphaRho, meanLT_subj,'MarkerFaceColor',color.general,'MarkerEdgeColor', color.general, 'MarkerFaceAlpha',M_alpha, 'MarkerEdgeAlpha',M_alpha); hold on
set(gca,'YLim', [0 25]) % small x axis, say in range of previous studies e.g. Garett & Daw. Check Lloyds studies too 
xlabel('Fit alpha')
ylabel('Patch leaving time (s)')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_ContrerasHuerta_fit_alpha_LT_correlation'],'-dsvg')

% statistics
[pearson_coeff, p_val]  = corr(minNLLFitParams.alphaRho, meanLT_subj)

%%%%% Kane 
load([dataPath 'subjectLT_kane'])
load([dataPath 'fitting_results_combined_M10_kane.mat'])

meanLT_subj = squeeze(mean(mean(subject_leave_times.mean(:,:,:),1),2));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
scatter(minNLLFitParams.alphaRho, meanLT_subj,'MarkerFaceColor',color.general,'MarkerEdgeColor', color.general, 'MarkerFaceAlpha',M_alpha, 'MarkerEdgeAlpha',M_alpha); hold on
set(gca,'YLim', [0 10]) % small x axis, say in range of previous studies e.g. Garett & Daw. Check Lloyds studies too 
xlabel('Fit alpha')
ylabel('Patch leaving time (s)')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis)
print([export_path 'fig5_Kane_fit_alpha_LT_correlation'],'-dsvg')

% statistics
[pearson_coeff, p_val]  = corr(minNLLFitParams.alphaRho, meanLT_subj)

