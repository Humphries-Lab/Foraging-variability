% model recovery of foraging models
% Emma Scholeyl, 17 March 2023
clearvars
close all;
addpath(genpath('~/Dropbox/foraging/code'))
model_table = readtable('~/Dropbox/foraging/details/input_model_table.xlsx');

model = [3 4 25 26];
Block = 1; % rich = 1, poor = 2
block_names = {'rich', 'poor'};
numParams  = 3;

% fmincon options
lowerBounds = [0,0,0];  % [alpha_Q, alpha_rho, beta] % parameter bounds
upperBounds = [1,1,100];  % arbitrary upper bound on beta to stop pathological behaviour
options = optimoptions('fmincon','Display','none'); % don't display

%% model recovery - plot confusion matrix
% extent to which simulated data from model A can be explained by model B

% set up the foraging environment
SetUpEnviron % get the environment parameters

CM = zeros(size(model,2));

for count = 1:100
    count

    figure(1); clf;
    FM = round(100*CM/sum(CM(1,:)))/100;
    t = imageTextMatrix(FM);
    set(t(FM'<0.3), 'color', 'w')
    hold on;
    [l1, l2] = addFacetLines(CM);
    set(t, 'fontsize', 22)
    title(['count = ' num2str(count)]);
    set(gca, 'xtick', [1:5], 'ytick', [1:5], 'fontsize', 28, ...
        'xaxislocation', 'top', 'tickdir', 'out')
    xlabel('fit model')
    ylabel('simulated model')

    drawnow

        testparams = [rand, rand, exprnd(3)]; % randomly select parameters each time

        % Model 3
        out = simulateM3_RLOn(testparams, Block, Env);
        SimData.Action = out.Action; 
        SimData.PatchOrder = Env.PatchOrder{Block};
        SimData.NextPatch = 0; % do this as the simulated data doesn't have stochastic patch prediction for this model, so need to tell models that do to start from scratch
        [BIC, iBEST, BEST] = fitAll(SimData, Env);
        CM(1,:) = CM(1,:) + BEST;

        % Model 4
        out = simulateM4_RLOn(testparams, Block, Env);
        SimData.Action = out.Action;
        SimData.PatchOrder = Env.PatchOrder{Block};
        SimData.NextPatch = 0;
        [BIC, iBEST, BEST] = fitAll(SimData, Env);
        CM(2,:) = CM(2,:) + BEST;

        % Model 25
        out = simulateM25_RLOn(testparams, Block, Env);
        SimData.Action = out.Action;
        SimData.PatchOrder = Env.PatchOrder{Block};
        SimData.NextPatch = out.NextPatch;
        [BIC, iBEST, BEST] = fitAll(SimData, Env);
        CM(3,:) = CM(3,:) + BEST;

        % Model 26
        out = simulateM26_RLOn(testparams, Block, Env);
        SimData.Action = out.Action;
        SimData.PatchOrder = Env.PatchOrder{Block};
        SimData.NextPatch = out.NextPatch;
        [BIC, iBEST, BEST] = fitAll(SimData, Env);
        CM(4,:) = CM(4,:) + BEST;

end
%
figure(1); 
title('')
set(gcf, 'Position', [811   417   500   400])
set(gca, 'fontsize', 28);


function [BIC, iBEST, BEST] = fitAll(Data, Env)

nStarts = 5; 
maxAlpha = [1 1];


% fit all of the models back to the simulated data of the model
[~, ~, BIC(1)] = fitM3_RLOn(Data,Env, nStarts, maxAlpha);
[~, ~, BIC(2)] = fitM4_RLOn(Data,Env, nStarts, maxAlpha);
[~, ~, BIC(3)] = fitM25_RLOn(Data,Env, nStarts, maxAlpha);
[~, ~, BIC(4)] = fitM26_RLOn(Data,Env, nStarts, maxAlpha);

% what's the best model for this data ?
[M, iBEST] = min(BIC);
BEST = BIC == M;
BEST = BEST / sum(BEST);

end


