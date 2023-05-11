% fitting RL model functions to foraging data
% Emma Scholey 9 Jun 2022
clear
close all;
addpath(genpath('~/Dropbox/foraging/code'))

model = 2;

SetUpEnviron % get the environment parameters
blockNames = {'rich', 'poor'};

options = optimoptions('fmincon','Display','none'); % don't display

% get their data
subj = [1:23 25:40];

% initialise
numSubjects = size(subj,2);
numBlocks = size(blockNames, 2);

load ~/Dropbox/foraging/raw_data/summary/young_variables/t_young.mat
group = {'Young_HC_ReDo/yc%g_forage.mat', subj};

%% fitting for each person in group with different starting points
nStarts = 50; % how many starting points to avoid local minima

names = {'rich', 'poor'};
numParams = 1;
minNLL = zeros([numSubjects, numBlocks]);
minNLLFitParams = zeros([numSubjects, numParams, numBlocks]);
FitParams = cell([numSubjects, 1]);
NLLEval = cell([numSubjects, 1]);

maxAlpha = [1 1]; % set max alpha to the range defined by best parameter fits (only need to recover parameters within this range)
% just set to maximum range (0-1)

timesteps = [0.25 0.5 1 2 3 6]; % timesteps to cover - related to travel time (6s)
for iT = 1:numel(timesteps)
    Env.TimeStep = timesteps(iT);
    iT
    for block = 1:numBlocks
        for k = 1:numSubjects
            SubjLT = t_young(t_young.subj == k,:); % extract their summarised leaving times
            Data = PrepSubjData(SubjLT, block, Env); % transform into state actions ready for fitting
            k
            [minNLL(k,block), minNLLFitParams(k,:,block)] = fitM2_MVT_RW(Data,Env, nStarts, maxAlpha);
        end
    end


    m(iT,:,block) = median(minNLLFitParams);

    figure; 
boxchart(squeeze(minNLLFitParams(:,1,:)))
title(sprintf('Best fit parameter values for beta, timestep = %d', timesteps(iT)), 'FontSize', 18)
xlabel('Environment', 'FontSize', 14); xticklabels({'Rich', 'Poor'});
ylabel('Beta value', 'FontSize', 14)
ylim([0 1.5])
end

 figure
 plot(timesteps, m(:,:,2))
 legend({'rich beta', 'poor beta'})


 %% manually for now 

 x = [0.125, 0.25, 0.5, 1, 2, 3]
 y = [0.26 0.21 0.18 0.13 0.1 0.07;0.37 0.3 0.25 0.19 0.14 0.105]

 figure 
 plot(x, y)
legend({'rich', 'poor'})
xlabel('Timestep (s)')
ylabel('Median best fit beta')
title('Fit betas decrease exponentially as timestep increases')
