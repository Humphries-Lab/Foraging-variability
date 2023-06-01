% model recovery of foraging models
% Emma Scholeyl, 17 March 2023
clearvars
close all;

allModels = 1:5;
modelNames = {'avgRR RW', 'patchRR RW', 'full RW', 'softmax only', 'patchRR RW fix beta'};

block = 1; % rich = 1, poor = 2
block_names = {'rich', 'poor'};
numParams  = 3;
numModels = size(allModels, 2);
% fmincon options
lowerBounds = [0,0,0];  % [alpha_Q, alpha_rho, beta] % parameter bounds
upperBounds = [1,1,100];  % arbitrary upper bound on beta to stop pathological behaviour
options = optimoptions('fmincon','Display','none'); % don't display

%% model recovery - plot confusion matrix
% extent to which simulated data from model A can be explained by model B

% set up the foraging environment
SetUpEnviron % get the environment parameters

CM = zeros(size(allModels,2));

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

    testparams = [rand, rand, exprnd(1)]; % randomly select parameters each time

    for m = 1:numModels
        sim_f = buildForagingModel(m, Env, [], block); % get function handle for this model
        out = sim_f(testparams); % simulate with the parameters
        out.PatchOrder = Env.PatchOrder{block}; % give unlimited patch order so that fitting back doesn't break 

        [BIC, iBEST, BEST] = fitAll(allModels, out, Env);
        CM(m,:) = CM(m,:) + BEST;

    end
    %
    figure(1);
    title('')
    set(gcf, 'Position', [811   417   500   400])
    set(gca, 'fontsize', 28);

end






function [BIC, iBEST, BEST] = fitAll(models, Data, Env)

numModels = size(models,2);

nStarts = 5;
options = optimoptions('fmincon','Display','none'); % don't display

for m = 1:numModels

    [~, NLL_f, paramsIndex] = buildForagingModel(m, Env, Data);

    lowerBounds = [0 0 0]; lowerBounds = lowerBounds(paramsIndex);
    upperBounds = [1 1 100]; upperBounds = upperBounds(paramsIndex);
    NLLEval = zeros([nStarts, 1]);


    % recover data
    parfor ii = 1:nStarts
        params0 = [defineBoundedParam(0, 1), defineBoundedParam(0, 1), exprnd(1)];
        params0 = params0(paramsIndex);
        [~,NLLEval(ii)] = fmincon(NLL_f,params0,[],[],[],[],lowerBounds, upperBounds,[],options);
    end
    numParams = sum(paramsIndex);
    numObservations = sum(Data.Action==2);
    minNLL = min(NLLEval);   % minimum negative log likelihood over all starting positions
    BIC(m) = numParams * log(numObservations) + 2*minNLL;

end

% what's the best model for this data ?
[M, iBEST] = min(BIC);
BEST = BIC == M;
BEST = BEST / sum(BEST);

end


