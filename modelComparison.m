%% model comparison of foraging
% Emma Scholey, 17 March 2023
clear  

%% user options 
study = 'leheron'; % study to simulate/fit data to

model = [1, 4:9]; % model numbers to compare 
modelTable = readtable('/Users/exs165/Dropbox/foraging/code/foragingModelTable.xlsx'); % change directory

blockPresentation = 'separate'; %% either 'combined' (fit as one continuous task) or 'separate' (fit rich and poor as separate blocks)

%% AIC/BIC

if contains(blockPresentation, 'combined')
    nBlock = 1; 
elseif contains(blockPresentation, 'separate')
    nBlock = 2; 
end 

nModels = size(model,2);

models_AIC = zeros(nModels, nBlock);
models_BIC = zeros(nModels, nBlock);

for m = 1:nModels
    load(sprintf('~/Dropbox/foraging/outputs/fitting/fitting_results_%s_M%d_new', blockPresentation, model(m)));
    ppts_AIC(:,:,m) = AIC;
    ppts_BIC(:,:,m) = BIC;
    models_AIC(m,:) = sum(AIC);
    models_BIC(m,:) = sum(BIC);
end

% compute posterior probabilities 
meanBICs = mean(models_BIC,2);
posteriorProbabilities = BICposterior(meanBICs);

% plot

figure
bar(meanBICs)
xticklabels(model)
ylabel('sum BIC')
title('BIC comparison of models (rich and poor)')


% find best model for each participant 
meanBIC = squeeze(mean(ppts_BIC,2));

bestModel = sum(meanBIC == min(meanBIC, [], 2))

bestModelRich = sum(squeeze(ppts_BIC(:,1,:)) == min(squeeze(ppts_BIC(:,1,:)), [], 2));

bestModelPoor = sum(squeeze(ppts_BIC(:,2,:)) == min(squeeze(ppts_BIC(:,2,:)), [], 2));