%% fitting models to foraging data
% Emma Scholey 9 Jun 2022
% latest update 30th May 2023

clear
close all;
addpath(genpath('~/Dropbox/foraging/code'))

%% user options

model = 4; %% see buildForagingModel.m for model list 
blockFlag = 'separate'; %%  either 'combined' (fit as one continuous task) or 'separate' (fit rich and poor as separate blocks) 

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

paramNames = {'alpha patch RR/Q', 'alpha rho', 'beta', 'lambda', 'epsilon', 'omega'};
paramNames = paramNames(paramsIndex); 

%% load and prepare subject data

% get subject data
subj = [1:23 25:40];
numSubjects = size(subj,2);

load ~/Dropbox/foraging/raw_data/summary/young_variables/t_young.mat
load ~/Dropbox/foraging/raw_data/experiencedAvgRR.mat

Data = cell([numSubjects, 1]);

for k = 1:numSubjects % for each subject
    for b = 1:numBlocks % for each block
        SubjLT = t_young(t_young.subj == k,:); % extract their summarised leaving times
        Data{k}{b} = PrepSubjData(SubjLT, b, Env, blockFlag); % transform into state actions ready for fitting
         Data{k}{b}.experiencedAvgRR{1} = experiencedAvgRR(k,1);
         Data{k}{b}.experiencedAvgRR{2} = experiencedAvgRR(k,2);
    end
end

% experiencedAvgRR = zeros([numSubjects, 2]);
% % log the experienced avgRR for each participant 
% for k = 1:numSubjects
%     for b = 1:numBlocks
%         experiencedAvgRR(k,1) = Data{k}{b}.experiencedAvgRR{1};
%         experiencedAvgRR(k,2) = Data{k}{b}.experiencedAvgRR{2};
%     end
% end
% save('~/Dropbox/foraging/raw_data/experiencedAvgRR.mat', 'experiencedAvgRR')

% assume no knowledge of environment for fitting
startValues = []; 

%% test NLL function on one participant
subjData = Data{8}{1}; % get example subject data
[NLL, out] = NLL_MVT_softmax([0.17], Env, subjData)

%% fitting for each person in group with different starting points

% fmincon options
nStarts = 50; % how many starting points to avoid local minima
lowerBounds = [0,0,0,-100, 0, -20];lowerBounds = lowerBounds(paramsIndex);  % {'alpha patch RR/Q', 'alpha rho', 'beta', 'lambda', 'epsilon', 'omega'} % parameter bounds
upperBounds = [1,1,10,300, 1, 20];upperBounds = upperBounds(paramsIndex);   % arbitrary upper bound on beta to stop pathological behaviour
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

        [~, NLL_f] = buildForagingModel(model, Env, subjData, startValues); % need to update NLL function with new subject data each time

         % Run fmincon
         parfor ii = 1:nStarts
             params0 = [rand, rand, exprnd(3), defineBoundedParam(0,50), rand, rand]; % choose full set of start parameters [alpha alpha beta lambda]
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
save_name = sprintf('~/Dropbox/foraging/outputs/fitting/fitting_results_%s_M%d', blockFlag, model);
save(save_name, 'AIC', 'BIC', 'minNLL', 'minNLLFitParams')

%% plots

close all
figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');

for i= 1:numParams
    nexttile;
    scatter(minNLLFitParams(:,i),minNLLFitParams(:,i+1))
    title(paramNames{i})
    xlabel(sprintf('Fit %s', paramNames{i}))
    xticklabels({'rich', 'poor'})
end
title(tl, 'Best fit parameter distributions')
