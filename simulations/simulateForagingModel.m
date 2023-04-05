%%simulateForagingModel

%% Set up
clear
close all
addpath(genpath('~/Dropbox/foraging/code'))
model = 4; % model type - see model table to check number to choose

NSim = 39; % how many runs
NReps = 10; % each participant simulated 10 times to reduce variability
% environment parameters
SetUpEnviron % generic script to set up foraging environment according to LeHeron et al (2020)

% simulate fresh
%load(sprintf('~/Dropbox/foraging/outputs/M%d/fitting_results_separate', model), 'minNLLFitParams');
%m = median(minNLLFitParams);
 %richParams = m(:,:,1);
 %poorParams = m(:,:,2);
% richParams = [0.7 0.1 4];
% poorParams = [0.7 0.1 4];
% blockParams = cat(3, repmat(richParams, [NSim, 1]), repmat(poorParams, [NSim, 1])); %[parameters for rich block; parameters for poor block % [alpha_Q, alpha_rho, beta]

% simulate after model fitting - using best fit parameters 
load(sprintf('~/Dropbox/foraging/outputs/M%d/fitting_results_separate', model), 'minNLLFitParams');
blockParams = repmat(minNLLFitParams, [NReps 1]); 

SimData = cell(NSim,1);
%% Run
for n = 1:NSim*NReps
    n
    for Block = 1:2 % for each environment: [rich, poor]
        params = blockParams(n,:,Block); 
        %out = simulateM1_MVT_RWRho(params, Block, Env); % change this depending on model to test
        %out = simulateM2_MVT_RW(params, Block, Env); % change this depending on model to test
        %out = simulateM3_RLOn(params, Block, Env); % change this depending on model to test
        out = simulateM4_RLOn(params, Block, Env); % change this depending on model to test
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
LTSEM = squeeze(std(LT))./sqrt(NSim*NReps);
LRRSEM = squeeze(std(LRR))./sqrt(NSim*NReps);

%% Store total reward stats
RewardMean = squeeze(mean(TotalReward))';

save(sprintf('~/Dropbox/foraging/outputs/M%d/SimulationResults.mat', model), 'RewardMean', "LTMean", "LTSEM", "LRRMean", "LRRSEM", "out", "NSim")

plotLeavingTimes(model)

%plotRun(model, SimData{1}{1}) % plot reward rates over the block - rich
%title('rich environment')
%plotRun(model, SimData{1}{2}) % poor
%title('poor environment')

%plotQ(model) % For RL Models only 