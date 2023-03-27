%%simulateForagingModel

%% Set up
clear
close all
addpath(genpath('~/Dropbox/foraging/code'))
model = 2; % model type - see model table to check number to choose

NSim = 100; % how many runs

% environment parameters
SetUpEnviron % generic script to set up foraging environment according to LeHeron et al (2020)

% simulate fresh
 richParams = [1    0.0004    9];
 poorParams = [1   0.0004    12];
 blockParams = cat(3, repmat(richParams, [NSim, 1]), repmat(poorParams, [NSim, 1])); %[parameters for rich block; parameters for poor block % [alpha_Q, alpha_rho, beta]

% simulate after model fitting - using best fit parameters 
%load(sprintf('~/Dropbox/foraging/outputs/M%d/fitting_results', model), 'minNLLFitParams');
%blockParams = minNLLFitParams; 

SimData = cell(NSim,1);
%% Run
for n = 1:NSim
    n
    for Block = 1:2 % for each environment: [rich, poor]
        params = blockParams(n,:,Block); 
        %out = simulateM1_MVT_RWRho(params, Block, Env); % change this depending on model to test
        out = simulateM2_MVT_RW(params, Block, Env); % change this depending on model to test
        %out = simulateM3_RLOn(params, Block, Env); % change this depending on model to test
        %out = simulateM4_RLOn(params, Block, Env); % change this depending on model to test
        %out = simulateM5_RLOn(params, Block, Env); % change this depending on model to test
        %out = simulateM21_RLOff(params, Block, Env); % change this depending on model to test
        %out = simulateM25_RLOn(params, Block, Env); % change this depending on model to test
        %out = simulateM26_RLOn(params, Block, Env); % change this depending on model to test

         SimData{n}{Block} = out;
         [LT(n,:,Block), LRR(n,:,Block), TotalReward(n,:,Block)] = summariseForaging(out); % get averages for each NSim
    end
end

%% compute means and standard errors across NSim
LTMean = squeeze(mean(LT));
LRRMean = squeeze(mean(LRR));
LTSEM = squeeze(std(LT))./sqrt(NSim);
LRRSEM = squeeze(std(LRR))./sqrt(NSim);

%% Store total reward stats
RewardMean = squeeze(mean(TotalReward))';

save(sprintf('~/Dropbox/foraging/outputs/M%d/SimulationResults.mat', model), 'RewardMean', "LTMean", "LTSEM", "LRRMean", "LRRSEM", "out", "NSim")

plotLeavingTimes(model)

plotRun(model, SimData{1}{1}) % plot reward rates over the block - rich
title('rich environment')
plotRun(model, SimData{1}{2}) % poor
title('poor environment')

%plotQ(model) % For RL Models only 