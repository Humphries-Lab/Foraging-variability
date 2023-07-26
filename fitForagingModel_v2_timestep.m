%% fitting models to foraging data - varying timestep
% Emma Scholey 13 July 2023

clear
close all;
addpath(genpath('~/Dropbox/foraging/code'))

%% user options

model = 11; %% see buildForagingModel.m for model list 
blockFlag = 'combined'; %%  either 'combined' (fit as one continuous task) or 'separate' (fit rich and poor as separate blocks) 

%% set up model and task 

if contains(blockFlag, 'combined')
    numBlocks = 1;
    blockNames = {'combined'};
elseif contains(blockFlag, 'separate')
    numBlocks = 2;
    blockNames = {'rich', 'poor'};
end

SetUpEnviron % get the environment parameters

[~, ~, paramsIndex] = buildForagingModel(model); % get the right parameters for the model 
numParams = sum(paramsIndex);

paramNames = {'alpha patch RR/Q', 'alpha rho', 'beta', 'lambda', 'epsilon', 'omega', 'timestep'};
paramNames = paramNames(paramsIndex); 

%% load and prepare subject data

% get subject data
subj = [1:23 25:40];
numSubjects = size(subj,2);

% assume no knowledge of environment for fitting
startValues = []; 

load ~/Dropbox/foraging/raw_data/summary/young_variables/t_young.mat
load ~/Dropbox/foraging/raw_data/experiencedAvgRR.mat

%% test NLL function on one participant
subjData.summary = t_young(t_young.subj == 8,:); % get example subject data
subjData.experiencedAvgRR = experiencedAvgRR(8,:);
[NLL, out] = NLL_MVT_softmax_dynamicbeta_timestep([6.7 1], Env, subjData);

%% fitting for each person in group with different starting points

% fmincon options
nStarts = 50; % how many starting points to avoid local minima
lowerBounds = [0,0,0,-100, 0, -20,0.01];lowerBounds = lowerBounds(paramsIndex);  % {'alpha patch RR/Q', 'alpha rho', 'beta', 'lambda', 'epsilon', 'omega'} % parameter bounds
upperBounds = [1,1,10,300, 1, 20,1];upperBounds = upperBounds(paramsIndex);   % arbitrary upper bound on beta to stop pathological behaviour
options = optimoptions('fmincon','Display','none'); % don't display

minNLL = zeros([numSubjects, numBlocks]);
minNLLFitParams = zeros([numSubjects, numParams, numBlocks]);
BIC = zeros([numSubjects, numBlocks]);
AIC = zeros([numSubjects, numBlocks]);

for block = 1:numBlocks

    for iS = 1:numSubjects

        iS
        SubjData.summary = t_young(t_young.subj == iS,:); % extract their summarised leaving times
        SubjData.experiencedAvgRR = experiencedAvgRR(iS,:);

        NLLEval = zeros([nStarts, 1]);
        FitParams = zeros([nStarts, numParams]);

        [~, NLL_f] = buildForagingModel(model, Env, SubjData, startValues); % need to update NLL function with new subject data each time

         % Run fmincon
         parfor ii = 1:nStarts
             params0 = [rand, rand, exprnd(3), defineBoundedParam(0,30), rand, rand, rand]; % choose full set of start parameters [alpha alpha beta lambda]
             params0 = params0(paramsIndex); % only take the parameters we need for this specific model

             [FitParams(ii,:),NLLEval(ii)] = fmincon(NLL_f,params0,[],[],[],[],lowerBounds, upperBounds,[],options);
         end

        % Find the best fitting parameter values
        minNLL = min(NLLEval);   % minimum negative log likelihood over all starting positions
        ix = find(minNLL == NLLEval);    % indices of location of minimum, to find the corresponding best fit parameters
        minNLLFitParams(iS,:,block) = FitParams(ix(1),:); % get corresponding parameter values at lowest NLL

        % number of observations for BIC/AIC depends on timestep used 
        [~, out] = NLL_MVT_softmax_dynamicbeta_timestep(minNLLFitParams(iS,:,block), Env, SubjData);

        % Calculate BIC/AIC
        BIC(iS,block) = numParams * log(out.numObservations) + 2*minNLL;
        AIC(iS,block) = 2/out.numObservations * minNLL + 2 * numParams/out.numObservations;

    end
end
clear FitParams NLLEval iS ix

m = median(minNLLFitParams);
save_name = sprintf('~/Dropbox/foraging/outputs/fitting/fitting_results_%s_M%d', blockFlag, model);
save(save_name, 'AIC', 'BIC', 'minNLL', 'minNLLFitParams')

%% plots

close all
figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');

for i= 1:numParams
    nexttile;
    boxchart(squeeze(minNLLFitParams(:,i,:)))
    title(paramNames{i})
    xlabel(sprintf('Fit %s', paramNames{i}))
    xticklabels({'rich', 'poor'})
end
title(tl, 'Best fit parameter distributions')
