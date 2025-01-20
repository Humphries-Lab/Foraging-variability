%% model comparison of foraging
% Emma Scholey, 17 March 2023

clear  
close all
addpath('./helperFunctions')

%% user options 
study = 'contrerashuerta'; % Options are leheron (field-human), contrerashuerta (berry-human), or kane (rat)

model = 1:5; % initial models
%model = [7 10 13 16]; % beta models
%model = [8 11 14 17]; % RW models 

modelTable = readtable('./foragingModelTable.xlsx'); 

%% AIC/BIC

nModels = size(model,2);

models_AIC = zeros(nModels, 1);
models_BIC = zeros(nModels, 1);

for m = 1:nModels
    load(sprintf('../data/fitting_data/fitting_results_M%d_%s', model(m),study), '-mat', 'minNLLFitParams','BIC', 'AIC');
    ppts_AIC(:,:,m) = AIC;
    ppts_BIC(:,:,m) = BIC;
    models_AIC(m,:) = sum(AIC);
    models_BIC(m,:) = sum(BIC);
end

% % compute posterior probabilities 
% posteriorProbabilities = BICposterior(squeeze(ppts_BIC));

% plot

figure
bar(models_BIC)
%xticklabels(model)
ylabel('sum BIC')
title('BIC comparison of models (rich and poor)')
ylim([min(models_BIC)-200,max(models_BIC)+200])

set(findall(gcf,'-property','FontSize'),'FontSize',18)

% find best model for each participant 
meanBIC = squeeze(mean(ppts_BIC,2));

subjectBestModelBIC = meanBIC == min(meanBIC, [], 2);
bestModel = sum(subjectBestModelBIC)
