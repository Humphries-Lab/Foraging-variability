%% model comparison of foraging
% Emma Scholey, 17 March 2023
clear  

%% user options 
study = 'leheron'; % study to simulate/fit data to

%model = [8, 9, 11, 12, 14, 15, 17, 18, 20, 21];
model = [1,26:29]; % model numbers to compare 
%model = [8:10]
modelTable = readtable('./foragingModelTable.xlsx'); 


%% AIC/BIC

nModels = size(model,2);

models_AIC = zeros(nModels, 1);
models_BIC = zeros(nModels, 1);

for m = 1:nModels
    load(sprintf('../data/fitting_data_240117/fitting_results_M%d_%s', model(m),study), '-mat', 'BIC', 'AIC');
    ppts_AIC(:,:,m) = AIC;
    ppts_BIC(:,:,m) = BIC;
    models_AIC(m,:) = sum(AIC);
    models_BIC(m,:) = sum(BIC);
end

% compute posterior probabilities 
sumBICs = sum(models_BIC,2);
%posteriorProbabilities = BICposterior(sumBICs);

sumAICs = sum(models_AIC,2);
% plot

figure
bar(sumBICs)
xticklabels(model)
ylabel('sum BIC')
title('BIC comparison of models (rich and poor)')
ylim([min(sumBICs)-200,max(sumBICs)+200])
xticklabels({'softmax (no c)', 'vary \beta vary c', 'vary \beta fix c', 'fix \beta vary c', 'fix \beta fix c'})

% find best model for each participant 
meanBIC = squeeze(mean(ppts_BIC,2));

subjectBestModelBIC = meanBIC == min(meanBIC, [], 2);
bestModel = sum(subjectBestModelBIC)
varNames = {'combined', 'modelN'};

%% formatting
fontsize = 14;
fontname = 'Arial';
htext = findobj(gca, 'type', 'text');
set(htext,'FontName',fontname,'FontSize',fontsize);


hAxes = findobj(gca, 'type', 'axes');
for i=1:numel(hAxes)
    set(get(hAxes(i),'XLabel'),'FontName',fontname,'FontSize',fontsize);
    set(get(hAxes(i),'YLabel'),'FontName',fontname,'FontSize',fontsize);
    set(hAxes(i),'FontName',fontname,'FontSize',fontsize);
end

% save for paper figures 
% models_BIC = array2table(models_BIC); models_BIC.modelN = model'; models_BIC.Properties.VariableNames = varNames;
% save_name = ['BIC_', blockPresentation, '_', study, '_',num2str(model), '.mat'];
% save_path = '../data/fig_data/';
% save([save_path, save_name],'models_BIC');