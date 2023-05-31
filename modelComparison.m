%% model comparison of foraging
% Emma Scholeyl, 17 March 2023


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

figure
t = tiledlayout('flow', 'TileSpacing', 'Compact');
title(t,'AIC comparison of models for each participant');
ylabel(t, 'AIC')

for subj = 1:size(ppts_AIC, 1)
    nexttile
    bar(squeeze(ppts_AIC(subj,:,:)));
end

sum(squeeze(ppts_AIC)'==min(squeeze(ppts_AIC)'),2)