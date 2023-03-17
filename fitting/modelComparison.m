%% model comparison of foraging
% Emma Scholeyl, 17 March 2023
clearvars
close all;
addpath(genpath('~/Dropbox/foraging/code'))
model_table = readtable('~/Dropbox/foraging/details/input_model_table.xlsx');

model = [3 4 25 26];
Block = 1; % rich = 1, poor = 2
block_names = {'rich', 'poor'};
numParams  = 3;

% fmincon options
lowerBounds = [0,0,0];  % [alpha_Q, alpha_rho, beta] % parameter bounds
upperBounds = [1,1,100];  % arbitrary upper bound on beta to stop pathological behaviour
options = optimoptions('fmincon','Display','none'); % don't display


%% AIC/BIC
models_AIC = zeros(size(model, 2), 2);
models_BIC = zeros(size(model, 2), 2);

for m = 1:size(model, 2)
    load(sprintf('~/Dropbox/foraging/outputs/M%d/fitting_results', model(m)));
    ppts_AIC(:,:,m) = AIC;
    ppts_BIC(:,:,m) = BIC;
    models_AIC(m,:) = sum(AIC);
    models_BIC(m,:) = sum(BIC);
end

figure
hist(mean(ppts_BIC(:,:,3),2))
% plot

figure
bar(mean(models_BIC,2))
xticklabels(model)
ylabel('BIC')
title('BIC comparison of models (rich and poor)')

figure
bar(mean(models_AIC,2))
xticklabels(model)
ylabel('AIC')
title('AIC comparison of models (rich and poor)')