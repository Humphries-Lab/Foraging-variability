%% fitting models to foraging data
% Emma Scholey 9 Jun 2022
% latest update 30th May 2023

clear
close all;
addpath(genpath('~/Dropbox/foraging/code'))

%% initialise

model = 5; %% CHANGE
blockFlag = 'combined'; %% CHANGE - either 'combined' (fit as one continuous task) or 'separate' (fit rich and poor as separate blocks) 

if contains(blockFlag, 'combined')
    numBlocks = 1;
elseif contains(blockFlag, 'separate')
    numBlocks = 2;
end

%% set up model and task 
SetUpEnviron % get the environment parameters
blockNames = {'rich', 'poor'};

[~, ~, paramsIndex] = buildForagingModel(model); % get the right parameters for the model 
numParams = sum(paramsIndex);

paramNames = {'alpha patch RR/Q', 'alpha rho', 'beta'};
paramNames = paramNames(paramsIndex); 

%% load and prepare subject data

% get subject data
subj = [1:23 25:40];
numSubjects = size(subj,2);

load ~/Dropbox/foraging/raw_data/summary/young_variables/t_young.mat

Data = cell([numSubjects, 1]);

for k = 1:numSubjects % for each subject
    for b = 1:numBlocks % for each block
        SubjLT = t_young(t_young.subj == k,:); % extract their summarised leaving times
        Data{k}{b} = PrepSubjData(SubjLT, b, Env, blockFlag); % transform into state actions ready for fitting
    end
end

%% fitting for each person in group with different starting points

% fmincon options
nStarts = 50; % how many starting points to avoid local minima
lowerBounds = [0,0,0];lowerBounds = lowerBounds(paramsIndex);  % [alpha_Q, alpha_rho, beta] % parameter bounds
upperBounds = [1,1,100];upperBounds = upperBounds(paramsIndex);   % arbitrary upper bound on beta to stop pathological behaviour
options = optimoptions('fmincon','Display','none'); % don't display

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

        [~, NLL_f] = buildForagingModel(model, Env, subjData); % need to update NLL function with new subject data each time

         % Run fmincon
         parfor ii = 1:nStarts
             params0 = [rand, rand, exprnd(3)]; % choose full set of start parameters [alpha alpha beta]
             params0 = params0(paramsIndex); % only take the parameters we need for this specific model

             [FitParams(ii,:),NLLEval(ii)] = fmincon(NLL_f,params0,[],[],[],[],lowerBounds, upperBounds,[],options);
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

