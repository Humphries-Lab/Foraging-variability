%% ---- Script to simulate Le Heron foraging task ---- %%
% Emma Scholey
% latest update 30th May 2023

clear
close all
addpath(genpath('~/Dropbox/foraging/code'))

%% user options
model = 10; % model type - see model table to check number to choose
blockFlag = 'combined'; %% either 'combined' (fit as one continuous task) or 'separate' (fit rich and poor as separate blocks)
simulateFitData = 1; % 1 if simulating subject fit data, 0 if want to simulate own parameters
NSim = 50; % this will override if simulating fit subject data 

%% set up agent and task
% environment parameters
blockNames = {'rich', 'poor'};
SetUpEnviron % generic script to set up foraging environment according to LeHeron et al (2020)
[~, ~, paramsIndex] = buildForagingModel(model);

if simulateFitData == 0
    if contains(blockFlag, 'separate') % parameters allowed to be different in each block
        richParams = [1, 0.67, 0.54, 3]; richParams = richParams(paramsIndex); % [alpha Q, alpha rho, beta, lambda]
        poorParams = [1, 0.0003, 0.22, 3]; poorParams = poorParams(paramsIndex);% [alpha Q, alpha rho, beta, lambda]

    elseif contains(blockFlag, 'combined')
        richParams = [0.5, 0.004, 0.5, 3]; richParams = richParams(paramsIndex);
        poorParams = richParams; % params have to be the same in both blocks
    end
    allParams = cat(3, repmat(richParams, [NSim, 1]), repmat(poorParams, [NSim, 1])); %[parameters for rich block; parameters for poor block % [alpha_Q, alpha_rho, beta]
    blockOrder = [repmat([2 1], [NSim/2, 1]); repmat([1 2], [NSim/2, 1])]; % do half [poor rich], half [rich poor]

elseif simulateFitData == 1
    load(sprintf('~/Dropbox/foraging/outputs/fitting/fitting_results_%s_M%d', blockFlag, model), 'minNLLFitParams');
    load('~/Dropbox/foraging/raw_data/BlockOrder.mat'); % load order of block presentation
    load('~/Dropbox/foraging/raw_data/experiencedAvgRR.mat'); % load their real average RR in the 2 blocks

    blockOrder = BlockOrder;
    NSim = 39;

    if contains(blockFlag, 'separate')
        allParams = minNLLFitParams;
    elseif contains(blockFlag, 'combined')
        allParams = cat(3, minNLLFitParams, minNLLFitParams); % duplicate for aligned block indexing later
    end

end

%% Run
SimData = cell(NSim,1);
LT = zeros(NSim, 3, 2);
LRR = zeros(NSim, 3, 2);

for n = 1:NSim
    n
   Env.experiencedAvgRR = experiencedAvgRR(n,:);

    for blockI = 1:2

        blockType = blockOrder(n,blockI); % which block do they see 1st or 2nd (to get order effects in simulations)
        Env.BlockPatchOrder = Env.PatchOrder{blockType}; % set the current patch order
        Env.CurrentBlock = blockType;
        
        if blockI == 2 && contains(blockFlag, 'combined') % we need to specify start values if it's the 2nd block and we're doing combined fitting of both blocks together
            startValues.Rho = out.Rho(end); % first value is last value from previous block
            startValues.EstimatedPatchRR = out.EstimatedPatchRR(end);
            startValues.sumReward = out.sumReward; 
            startValues.timeInPreviousBlock = Env.BlockTime;
        else
            startValues.Rho = 0;
            startValues.EstimatedPatchRR = 0;
            startValues.sumReward = 0;
            startValues.timeInPreviousBlock = 0; 
        end

        params = allParams(n,:,blockType);

        sim_f = buildForagingModel(model, Env, [], startValues); % get function handle for this model
        out = sim_f(params); % simulate with the parameters

        SimData{n}{blockType} = out;

        % extract leaving times for each patch type
        for ii = 1:3
            LT(n,ii,blockType) = mean(out.LeavingTime(out.PatchOrder == ii));
            LRR(n,ii,blockType) = mean(out.LeavingRR(out.PatchOrder == ii));
        end

        % how much reward did they get?
        TotalReward(n,blockType) = sum(out.Reward);
    end
end

%% compute means and standard errors across NSim
LTMean = squeeze(mean(LT));
LRRMean = squeeze(mean(LRR));
LTSEM = squeeze(std(LT))./sqrt(NSim);
LRRSEM = squeeze(std(LRR))./sqrt(NSim);
RewardMean = squeeze(mean(TotalReward))';

%save(sprintf('~/Dropbox/foraging/sim_data/simData_M%d.mat', model), 'RewardMean', "LTMean", "LTSEM", "out", "NSim", "SimData")

%% figure - plot mean leaving times, compared to real subject data
close all

load('~/Dropbox/foraging/raw_data/summary/young_variables/subMean_young.mat') % load real data

% set up colours
c = [0.6353    0.0784    0.1843; 0    0.4471    0.7412]; % [rich, poor]

% plot the actual experimental data
NSub = size(subMean_young, 3);
SubLT(:,:,1) = squeeze(subMean_young(1,:,:)); % rich environment
SubLT(:,:,2) = squeeze(subMean_young(2,:,:)); % poor environment
SubLTMean = [mean(SubLT(:,:,1), 2), mean(SubLT(:,:,2), 2)];
SubSEM = [std(SubLT(:,:,1),[],2)./sqrt(NSub), std(SubLT(:,:,2), [], 2)./sqrt(NSub)]; % Standard Error

ts = tinv([0.025  0.975],NSub-1);      % T-Score
SubCI = ts(2).*SubSEM; % CIs

figure()
hold on
title('Mean patch leaving times');
xlabel('Patch yield type')
ylabel('Time in patch (s)')

errorbar([1 2 3],SubLTMean(:,1), SubCI(:,1),'LineStyle','--', 'Color',c(1,:) ,'LineWidth', 1);
hold on
errorbar([1 2 3],SubLTMean(:,2), SubCI(:,2),'LineStyle','--', 'Color',c(2,:) ,'LineWidth', 1);

%overlay simulated data
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

figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');

for i= 1:2
    ax = nexttile;
    plot(1:Env.BlockTime, SimData{8}{i}.PatchRR, 'LineWidth', 1) % plot patch reward rate
    hold on
    plot(1:Env.BlockTime, SimData{8}{i}.Rho, 'LineWidth', 1) % plot estimated average RR
    legend('Patch RR', 'Experienced average RR', 'FontSize', 16, 'FontName', 'Helvetica')
    xlabel('Time (s)','FontSize', 16, 'FontName', 'Helvetica');
    ylabel('Reward rates','FontSize', 16, 'FontName', 'Helvetica');
    title(sprintf('%s environment', blockNames{i}));
end

% plot how betas change over course of experiment (for dynamic beta models
% only) % for each parameter, plot sim vs fit
% figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');
% for i= 1:2
%     ax = nexttile;
%     plot(1:Env.BlockTime, SimData{1}{i}.Beta, 'LineWidth', 1) % plot patch reward rate
%     %hold on
%     %plot(1:Env.BlockTime, SimData{2}{i}.experiencedAvgRR, 'LineWidth', 1) % plot estimated average RR
%     legend('Beta', 'Experienced average RR', 'FontSize', 16, 'FontName', 'Helvetica')
%     xlabel('Time (s)','FontSize', 16, 'FontName', 'Helvetica');
%     title(sprintf('%s environment', blockNames{i}));
%     %ylim([0,1])
% end

%% compare simulated leaving times with real leaving times for each participant 

load ~/Dropbox/foraging/raw_data/summary/young_variables/t_young.mat

figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');

for iSubj = 1:NSub
    ax = nexttile;
    simLT = SimData{iSubj}{i}.LeavingTime;
    subjLT = t_young.leaveT(t_young.subj == iSubj & t_young.env == 1);

    plot(simLT, 'LineWidth', 1); hold on; plot(subjLT, 'LineWidth', 1)
    ylim([0 70])
end
