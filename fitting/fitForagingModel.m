% fitting RL model functions to foraging data
% Emma Scholey 9 Jun 2022
clear
close all;
addpath(genpath('~/Dropbox/foraging/code'))

model = 2;

SetUpEnviron % get the environment parameters
blockNames = {'rich', 'poor'};

% fmincon options
lowerBounds = [0,0,0];  % [alpha_Q, alpha_rho, beta] % parameter bounds
upperBounds = [1,1,100];  % arbitrary upper bound on beta to stop pathological behaviour
options = optimoptions('fmincon','Display','none'); % don't display

% get their data
subj = [1:23 25:40];

% initialise
numSubjects = size(subj,2);
numParams = size(lowerBounds,2);
numBlocks = size(blockNames, 2);

%subj = 1:4;
load ~/Dropbox/foraging/raw_data/summary/young_variables/t_young.mat
group = {'Young_HC_ReDo/yc%g_forage.mat', subj};

Data = cell([numSubjects, 1]);

for k = 1:numSubjects % for each subject
    for b = 1:numBlocks % for each block [rich poor]
        SubjLT = t_young(t_young.subj == k,:); % extract their summarised leaving times
        Data{k}{b} = PrepSubjData(SubjLT, b, Env); % transform into state actions ready for fitting
    end
end

clear group k b SubjLT t_young

%% grid search for subjects - marginal likelilood distributions
% range of values for grid search of starting points for search
AlphaVector = 0:0.05:1;
BetaVector = 0:2.5:50;

MinNLLEval = cell([numSubjects, 1]);
ParamsMinNLLEval = cell([numSubjects, 1]);
NLLEval = cell([numSubjects, 1]);

for k = 1:numSubjects
    k
    for block = 1:numBlocks % for each environment: [rich, poor]
        [MinNLLEval{k}{block}, ParamsMinNLLEval{k}{block}, NLLEval{k}{block}] = GridSearchForaging(Data{k}{block}, Env, AlphaVector, BetaVector); % change this depending on model to test
    end
end

% plot marginal likelihoods across grid for each subject
plotMarginalLikelihoods(NLLEval, AlphaVector,BetaVector)


%% fitting for each person in group with different starting points
nStarts = 50; % how many starting points to avoid local minima

names = {'rich', 'poor'};

minNLL = zeros([numSubjects, numBlocks]);
minNLLFitParams = zeros([numSubjects, numParams, numBlocks]);
%minNLLFitParamsSE = zeros([numSubjects, numParams, numBlocks]);
FitParams = cell([numSubjects, 1]);
NLLEval = cell([numSubjects, 1]);

BIC = zeros([numSubjects, numBlocks]);
AIC = zeros([numSubjects, numBlocks]);
maxAlpha = [upperBounds(1) upperBounds(2)]; % what range of starting values should we look over for Alphas (this setting is more important for model fitting)
% just set to maximum range (0-1)

for block = 1:numBlocks
    for k = 1:numSubjects
        k
        [minNLL(k,block), minNLLFitParams(k,:,block), BIC(k,block), AIC(k,block), FitParams{k}{block}, NLLEval{k}{block}] = fitM2_MVT_RW(Data{k}{block},Env, nStarts, maxAlpha);
        %[minNLL(k,block), minNLLFitParams(k,:,block), BIC(k,block), AIC(k,block), FitParams{k}{block}, NLLEval{k}{block}] = fitM5_RLOn(Data{k}{block},Env, nStarts, maxAlpha);
        %[minNLL(k,block), minNLLFitParams(k,:,block), BIC(k,block), AIC(k,block), FitParams{k}{block}, NLLEval{k}{block}] = fitM21_RLOff(Data{k}{block},Env, nStarts, maxAlpha);
        %[minNLL(k,block), minNLLFitParams(k,:,block), BIC(k,block), AIC(k,block), FitParams{k}{block}, NLLEval{k}{block}] = fitM25_RLOn(Data{k}{block},Env, nStarts, maxAlpha);
    end
    save_name = sprintf('~/Dropbox/foraging/outputs/M%d/fitting_results_%s', model, names{block});
    save(save_name, 'AIC', 'BIC', 'minNLL', 'minNLLFitParams', 'FitParams')
end

% plots

plotParamDist(FitParams) % plot distribution of parameter values for each subject (and each parameter, for each environment)

m = median(minNLLFitParams);
plotBestParams(minNLLFitParams) % plot best parameters across subjects for each block (and median)

%% parameter recovery

Block = 1;

nSim = 250; % number of iterations
maxAlpha = [0.85 0.5]; % set max alpha to the range defined by best parameter fits (only need to recover parameters within this range)

SimParams = zeros(numParams, nSim);
FitParams = zeros(numParams, nSim);

for ii = 1:nSim
    ii
    TestParams = [0.001+rand(1,1)*(maxAlpha(1)-0.001), 0.001+rand(1,1)*(maxAlpha(2)-0.001), exprnd(10)]; % randomly select parameters each time

    % simulate data
    out = simulateM2_MVT_RW(TestParams, Block, Env); % change this depending on model to test
    %out = simulateM5_RLOn(TestParams, Block, Env); % change this depending on model to test
    %out = simulateM21_RLOff(TestParams, Block, Env); % change this depending on model to test
    %out = simulateM25_RLOn(TestParams, Block, Env); % change this depending on model to test

    SimData.Action = out.Action;
    SimData.PatchOrder = Env.PatchOrder{Block};
    if model > 3 % for every model apart from MVT 
        SimData.NextPatch = out.NextPatch;
    end
    % recover data
    nStarts = 5; % how many starting locations to begin fitting (avoid local minima)
    [~, minNLLFitParams, ~, ~, SimFitParams] = fitM2_MVT_RW(SimData,Env, nStarts, maxAlpha);
    %[~, minNLLFitParams, ~, ~, SimFitParams] = fitM5_RLOn(SimData,Env, nStarts, maxAlpha);
    %[~, minNLLFitParams, ~, ~, SimFitParams] = fitM21_RLOff(SimData,Env, nStarts, maxAlpha);
    %[~, minNLLFitParams, ~, ~, SimFitParams] = fitM25_RLOn(SimData,Env, nStarts, maxAlpha);

    for iParam = 1:numParams
        SimParams(iParam,ii) = TestParams(iParam);
        FitParams(iParam,ii) = minNLLFitParams(iParam);
    end
end

plotParamRecovery(SimParams, FitParams) %  compute and plot correlations between simulated/fit parameters, and between fit parameters (trade-off)

%% Checks for parameter recovery - extreme states
for ii = 1:nSim
    max_Q_length(ii) = max(sum(Q_stay{ii} ~= 0));
end

residuals_alpha_patch = sqrt((FitParams(1,:) - SimParams(1,:)).^2);
residuals_alpha_rho = sqrt((FitParams(2,:) - SimParams(2,:)).^2);
residuals_beta = sqrt((FitParams(3,:) - SimParams(3,:)).^2);

% lower distance should correlate with shorter Q length
corrcoef(log(max_Q_length), log(residuals_alpha_patch))
corr(max_Q_length', residuals_alpha_patch', 'Type', 'Spearman')

corrcoef(log(max_Q_length), log(residuals_alpha_rho))
corr(max_Q_length', residuals_alpha_rho', 'Type', 'Spearman')

corrcoef(log(max_Q_length), log(residuals_beta))
corr(max_Q_length', residuals_beta', 'Type', 'Spearman')

figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');
ax = nexttile;
plot(log(residuals_alpha_patch), log(max_Q_length), 'o', 'markersize', 8, 'linewidth', 1)
title('alpha patch fit residuals vs max Q stay length (Extreme states)')
xlabel('alpha patch residuals')
ylabel('Length of extreme state')
set(ax,'xscale' ,'log')
set(ax,'yscale' ,'log')

ax = nexttile;
plot(log(residuals_alpha_rho), log(max_Q_length), 'o', 'markersize', 8, 'linewidth', 1)
title('alpha rho fit residuals vs max Q stay length (Extreme states)')
xlabel('alpha rho residuals')
ylabel('Length of extreme state')
set(ax,'xscale' ,'log')
set(ax,'yscale' ,'log')

ax = nexttile;
plot(log(residuals_beta), log(max_Q_length), 'o', 'markersize', 8, 'linewidth', 1)
title('beta fit residuals vs max Q stay length (Extreme states)')
xlabel('beta residuals')
ylabel('Length of extreme state')
set(ax,'xscale' ,'log')
set(ax,'yscale' ,'log')

%% checks - simulation and fitting scripts end up with the same outputs
Block = 1;

TestParams = [0.4, 0.1, 5];

% simulate data
out = simulateM2_MVT_RW(TestParams, Block, Env); % change this depending on model to test
%out = simulateM5_RLOn(TestParams, Block, Env); % change this depending on model to test
%out = simulateM21_RLOff(TestParams, Block, Env); % change this depending on model to test
%out = simulateM25_RLOn(TestParams, Block, Env); % change this depending on model to test
SimData.Action = out.Action;
SimData.PatchOrder = Env.PatchOrder{Block};
if model > 3 % for every model apart from MVT
    SimData.NextPatch = out.NextPatch;
end

[SubjNLLEval, RecoveredData] = NLL_M2_MVT_RW(TestParams, Env, SimData);
%[SubjNLLEval, RecoveredData] = NLL_M5_RLOn(TestParams, Env, SimData);
%[SubjNLLEval, RecoveredData] = NLL_M21_RLOff(TestParams, Env, SimData);
%[SubjNLLEval, RecoveredData] = NLL_M25_RLOn(TestParams, Env, SimData);

% issue is next patch - Q values start to diverge as soon as there's a
% difference in the NextPatch prediction. I've now changed the next patch
% prediction to be just mode (not stochastic). If this improves recovery,
% then will need to log simulated next patch

%% checks - distance from start parameter to fit parameter
% showing good convergence for most subjects, suggests fitting as expected.
%

% plot distance from start to fit parameters
nStarts = 50; % how many starting points to avoid local minima

maxAlpha = [upperBounds(1) upperBounds(2)]; % what range of starting values should we look over for Alphas (this setting is more important for model fitting)
% just set to maximum range (0-1)

for block = 1:numBlocks
    for k = 1:numSubjects
        k
        [~, ~, ~, ~, FitParams, ~, StartParams] = Fit_M25_RLOn(Data{k}{block},Env, nStarts, maxAlpha);
        plotFitDistance(FitParams, StartParams, nStarts)
        pause
    end
end

