%% fitting models to foraging data
% Emma Scholey 9 Jun 2022
% latest update 30th May 2023

clear
close all;
addpath(genpath('~/Dropbox/foraging/code'))

%% initialise

model = 25;

SetUpEnviron % get the environment parameters
blockNames = {'rich', 'poor'};

%% set up model parameter options for later
% ----- list of models ----- %
% M1 - MVT proper
% M2 - MVT RW   agent.alpha patch, agent.alpha rho, agent.beta
% M3 - RL on
% M4 - RL on
% M5 - RL on

% M21 RL off
% M25 RL on
% M26 RL on

% define which parameters we need for this model - index is model x
% parameter
modelParamsIndex = logical([1 1 1]); % M2
%modelParamsIndex = logical(paramIndex(model,:));
paramNames = {'alpha patch RR/Q', 'alpha rho', 'beta'};

%% load and prepare subject data

% get subject data
subj = [1:23 25:40];
numSubjects = size(subj,2);
numBlocks = size(blockNames, 2);

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

%% fitting for each person in group with different starting points

% fmincon options
nStarts = 50; % how many starting points to avoid local minima
lowerBounds = [0,0,0];  % [alpha_Q, alpha_rho, beta] % parameter bounds
upperBounds = [1,1,100];  % arbitrary upper bound on beta to stop pathological behaviour
options = optimoptions('fmincon','Display','none'); % don't display

numBlocks = size(blockNames, 2);

numParams = sum(modelParamsIndex);

minNLL = zeros([numSubjects, numBlocks]);
minNLLFitParams = zeros([numSubjects, numParams, numBlocks]);
BIC = zeros([numSubjects, numBlocks]);
AIC = zeros([numSubjects, numBlocks]);

for block = 1:numBlocks
    for iS = 1:numSubjects
        iS
        subjData = Data{iS}{block};
        numObservations = sum(subjData.Action == 2); % only count states when subject stays in patch (this is when log likelihood updates as they make a decision in that state) 

         NLLEval = zeros([nStarts, 1]);
         FitParams = zeros([nStarts, numParams]);
        % Run fmincon
        parfor ii = 1:nStarts
            params0 = [rand, rand, exprnd(3)]; % choose start parameters [alpha alpha beta]
            params0 = params0(modelParamsIndex); % only take the parameters we need for this model

            %f = @(x)NLL_M1_MVT(x, Env, subjData, block);  % need to tell MVT proper which block, to get optimal average RR as threshold
            %f = @(x)NLL_M2_MVT_RW(x, Env, subjData);
            %f = @(x)NLL_M3_RLOn(x, Env, subjData);
            %f = @(x)NLL_M4_RLOn(x, Env, subjData);
            %f = @(x)NLL_M5_RLOn(x, Env, subjData);

            %f = @(x)NLL_M21_RLOff(x, Env, subjData);
            f = @(x)NLL_M25_RLOn(x, Env, subjData);
            %f = @(x)NLL_M26_RLOn(x, Env, subjData);

            [FitParams(ii,:),NLLEval(ii)] = fmincon(f,params0,[],[],[],[],lowerBounds, upperBounds,[],options);
        end
        % Find the best fitting parameter values
        minNLL = min(NLLEval);   % minimum negative log likelihood over all starting positions
        ix = find(minNLL == NLLEval);    % indices of location of minimum, to find the corresponding best fit parameters
        minNLLFitParams(iS,:,block) = FitParams(ix(1),:); % get corresponding parameter values at lowest NLL

        % Calculate BIC/AIC
        BIC(iS,block) = numParams * log(numObservations) + 2*minNLL;
        AIC(iS,block) = 2/numObservations * minNLL + 2 * numParams/numObservations;

    end
end
clear FitParams NLLEval iS ix

m = median(minNLLFitParams);
save_name = sprintf('~/Dropbox/foraging/outputs/fitting/fitting_results_separate_M%d', model);
save(save_name, 'AIC', 'BIC', 'minNLL', 'minNLLFitParams')

%% plots

close all
figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');

for i= 1:numParams
    nexttile;
    boxchart(squeeze(minNLLFitParams(:,i,:)))
    title(paramNames{i})
    xlabel(sprintf('Fit %s', paramNames{i}))
    xticklabels(blockNames)
end
title(tl, 'Best fit parameter distributions')

