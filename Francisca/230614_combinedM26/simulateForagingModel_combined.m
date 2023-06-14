%% ---- Script to simulate Le Heron foraging task ---- %%
% Emma Scholey 
% latest update 30th May 2023 

clear
close all
%addpath(genpath('./foraging/code'))

%% initialise
model = 26; % model type - see model table to check number to choose

NSim = 100; % how many runs - set to participant number (39) or arbitrary value if simulating from scratch 

% environment parameters
SetUpEnviron_Francisca % generic script to set up foraging environment according to LeHeron et al (2020)

% agent parameters 

% simulate from scratch
params = [0.14 0.04 3.7]; % [alpha Q, alpha rho, beta]
simParams = repmat(params, [NSim, 1]); %[alpha_Q, alpha_rho, beta]
blockOrder = [repmat([2 1], [NSim/2, 1]); repmat([1 2], [NSim/2, 1])]; 

% OR simulate after model fitting - using best fit parameters 
%load(sprintf('~/Dropbox/foraging/outputs/fitting/fitting_results_separate_M%d', model), 'minNLLFitParams');
%blockParams = minNLLFitParams; 

SimData = cell(NSim,1);

%% Run
for n = 1:NSim
    n
    params = simParams(n,:);
    out = simulateM26_RLOn_combined(params, blockOrder(n,:), Env); % change this depending on model to test

    SimData{n} = out;

    % extract leaving times for each patch type
    for ii = 1:3
        for block = 1:2
            LT(n,ii,block) = mean(out.LeavingTime{block}(out.PatchOrder{block}(1:end-1) == ii), 'omitnan');
        end
    end
end


%% compute means and standard errors across NSim
LTMean = squeeze(mean(LT));
LTSEM = squeeze(std(LT))./sqrt(NSim);

%% figure - plot mean leaving times, compared to real subject data 
% close all
% 
% load('./foraging/raw_data/summary/young_variables/subMean_young.mat') % load real data
% 
% set up colours
c = [0.6353    0.0784    0.1843; 0    0.4471    0.7412]; % [rich, poor]
% 
% % plot the actual experimental data
% NSub = size(subMean_young, 3);
% SubLT(:,:,1) = squeeze(subMean_young(1,:,:)); % rich environment
% SubLT(:,:,2) = squeeze(subMean_young(2,:,:)); % poor environment
% SubLTMean = [mean(SubLT(:,:,1), 2), mean(SubLT(:,:,2), 2)];
% SubSEM = [std(SubLT(:,:,1),[],2)./sqrt(NSub), std(SubLT(:,:,2), [], 2)./sqrt(NSub)]; % Standard Error
% 
% ts = tinv([0.025  0.975],NSub-1);      % T-Score
% SubCI = ts(2).*SubSEM; % CIs 

figure()
hold on
title('Mean patch leaving times');
xlabel('Patch yield type')
ylabel('Time in patch (s)')

% errorbar([1 2 3],SubLTMean(:,1), SubCI(:,1),'LineStyle','--', 'Color',c(1,:) ,'LineWidth', 1);
% hold on
% errorbar([1 2 3],SubLTMean(:,2), SubCI(:,2),'LineStyle','--', 'Color',c(2,:) ,'LineWidth', 1);

% overlay simulated data
%load(sprintf('~/Dropbox/foraging/sim_data/simData_M%d.mat', model)) % load simulated results for each model

ts = tinv([0.025  0.975],NSim-1);      % T-Score
LTCI = LTSEM .* ts(2);

errorbar([1 2 3], LTMean(:,1), LTCI(:,1), 'LineWidth', 2, 'Color', c(1,:)); % rich environment
errorbar([1 2 3], LTMean(:,2), LTCI(:,2), 'LineWidth', 2, 'Color', c(2,:)); % poor environment


title(sprintf('M%d', model))

% tidy up axes
ax = gca;
xlim(ax, xlim(ax) + [-1,1]*range(xlim(ax)).* 0.1)
xticks(ax, [1 2 3]); xticklabels(ax, {'Low'; 'Medium'; 'High'}); xlabel('Patch quality')
ylim([0, 30]); ylabel('Patch leaving times (s)')

legend({'Subject - rich', 'Subject - poor', 'Agent - rich', 'Agent - poor'}, 'Location', 'north');
subtitle(sprintf('Error bars show 95%% CI. NSim = %d', NSim))

%% plot behaviour from a single run

figure
plot(1:Env.TaskTime+1, out.PatchRR, 'LineWidth', 1) % plot patch reward rate
hold on
plot(1:Env.TaskTime+1, out.Rho, 'LineWidth', 1) % plot estimated average RR
legend('Patch RR', '\rho: Estimated average RR', 'FontSize', 16, 'FontName', 'Helvetica')
xlabel('Time (s)','FontSize', 16, 'FontName', 'Helvetica');
ylabel('Reward rates','FontSize', 16, 'FontName', 'Helvetica');
