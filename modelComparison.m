%% model comparison of foraging
% Emma Scholey, 17 March 2023
clear  

%% user options 
model = 6; 
modelNames = {'M6'}; % automate this eventually
blockFlag = 'combined'; %% either 'combined' (fit as one continuous task) or 'separate' (fit rich and poor as separate blocks)



%% AIC/BIC

if contains(blockFlag, 'combined')
    nBlock = 1; 
elseif contains(blockFlag, 'separate')
    nBlock = 2; 
end 

nModels = size(model,2);

models_AIC = zeros(nModels, nBlock);
models_BIC = zeros(nModels, nBlock);

for m = 1:nModels
    load(sprintf('~/Dropbox/foraging/outputs/fitting/fitting_results_%s_M%d', blockFlag, model(m)));
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
xticklabels(modelNames)
ylabel('sum BIC')
title('BIC comparison of models (rich and poor)')

figure
bar(meanBICs)
xticklabels(modelNames)
ylabel('sum AIC')
title('AIC comparison of models (rich and poor)')


% find best model for each participant 
meanBIC = squeeze(mean(ppts_BIC,2));

bestModel = sum(meanBIC == min(meanBIC, [], 2))

bestModelRich = sum(squeeze(ppts_BIC(:,1,:)) == min(squeeze(ppts_BIC(:,1,:)), [], 2));

bestModelPoor = sum(squeeze(ppts_BIC(:,2,:)) == min(squeeze(ppts_BIC(:,2,:)), [], 2));


% figure
% t = tiledlayout('flow', 'TileSpacing', 'Compact');
% title(t,'AIC comparison of models for each participant');
% ylabel(t, 'AIC')
% 
% sum(squeeze(ppts_AIC)'==min(squeeze(ppts_AIC)'),2)