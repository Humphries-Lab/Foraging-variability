%% model comparison of foraging
% Emma Scholey, 17 March 2023
clear  
close all
addpath('./helperFunctions')

%% user options 
study = 'contrerashuerta'; % study to simulate/fit data to

%model = [8, 9, 11, 12, 14, 15, 17, 18, 20, 21];
model = [3,7]; % model numbers to compare 
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
posteriorProbabilities = BICposterior(squeeze(ppts_BIC));

% plot

figure
bar(models_BIC)
xticklabels(model)
ylabel('sum BIC')
title('BIC comparison of models (rich and poor)')
ylim([min(models_BIC)-200,max(models_BIC)+200])
%xticklabels({'vary \beta vary c', 'vary \beta fix c', 'fix \beta vary c', 'fix \beta fix c'})
%xticklabels({'softmax', 'vary \beta vary c', 'vary \beta fix c'})
set(findall(gcf,'-property','FontSize'),'FontSize',18)

% find best model for each participant 
meanBIC = squeeze(mean(ppts_BIC,2));

subjectBestModelBIC = meanBIC == min(meanBIC, [], 2);
bestModel = sum(subjectBestModelBIC)

% save for paper figures 
 models_BIC = array2table(models_BIC); models_BIC.modelN = model';
%  save_name = ['BIC_', study, '_',num2str(model), '.mat'];
%  save_path = '../data/fig_data/fig3/';
%  save([save_path, save_name],'models_BIC', 'posteriorProbabilities');