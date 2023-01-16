%% validateForagingModel
% clear
model = 25; % model type - see model table to check number to choose

load(sprintf('~/Dropbox/foraging-project/results/M%d/fitting_results_rich.mat', model)) % load best fit parameters from real data
params{1} = minNLL_params_fitted;
load(sprintf('~/Dropbox/foraging-project/results/M%d/fitting_results_poor.mat', model)) % load best fit parameters from real data
params{2} = minNLL_params_fitted;
NSim = size(minNLL_params_fitted,1); % simulate each person 

Data = cell(NSim,1);
%% Run
for n = 1:NSim
    n
    for Block = 1:2 % for each environment: [rich, poor]
        out = simulateM25_RLOn(params, Block); % change this depending on model to test
        Data{n}{Block} = out;
        [LT(n,:,Block), LRR(n,:,Block), TotalReward(n,:,Block)] = summariseForaging(out); % get averages for each NSim
    end
end

%% compute means across NSim
LTMean = squeeze(mean(LT));
LRRMean = squeeze(mean(LRR));
LTStandardDev= squeeze(std(LT));
LRRStandardDev = squeeze(std(LRR));

%% Store total reward stats
RewardMean = squeeze(mean(TotalReward))';

save(sprintf('~/Dropbox/foraging/outputs/M%d/ValidationResults.mat', model), 'RewardMean', "LTMean", "LTStandardDev", "LRRMean", "LRRStandardDev")

