% model recovery of foraging models
% Emma Scholeyl, 17 March 2023
clear
close all

addpath('./helperFunctions')
addpath('./analyticalComputation')
%% user options
allModels = 1:5;
blockFlag = 'combined'; %% either 'combined' (fit as one continuous task) or 'separate' (fit rich and poor as separate blocks)
modelNames = {'avgRR RW', 'patchRR RW', 'full RW', 'softmax only', 'patchRR RW fix beta'};
block = 1; % rich = 1, poor = 2. Only required for separate block model fitting

%% model recovery - plot confusion matrix
% extent to which simulated data from model A can be explained by model B
block_names = {'rich', 'poor'};
numModels = size(allModels, 2);
options = optimoptions('fmincon','Display','none'); % don't display

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

    for m = 1:numModels

        %% simulate data
        simData = table; % container table for simulated data
        [~, ~, paramsIndex] = buildForagingModel(m); % find which parameters this model uses
        params = [defineBoundedParam(0, 1), defineBoundedParam(0, 0.05), exprnd(1)]; % random parameters each time
        params = params(paramsIndex); % restrict to only the parameters the models needs

        for blockI = 1:2

            blockType = blockI; % for recovery, assume always seen rich block first (same order as blockI)
            Env.BlockPatchOrder = Env.PatchOrder{blockType}; % set the current patch order

            if blockI == 2 && contains(blockFlag, 'combined') % we need to specify start values if it's the 2nd block and we're doing combined fitting of both blocks together
                startValues.Rho = out.Rho(end); % first value is last value from previous block
                startValues.EstimatedPatchRR = out.EstimatedPatchRR(end);
            else
                startValues.Rho = 0;
                startValues.EstimatedPatchRR = 0;
            end

            sim_f = buildForagingModel(m, Env, [], startValues); % get function handle for this model
            out = sim_f(params); % simulate with the parameters
            % make a table of simulated data, similar to the real participant
            % data (whilst it's a bit clunky, it helps for understanding if/how NLL is working differently to
            % simulation script)
            tmp = table;
            tmp.leaveT = out.LeavingTime;
            tmp.patch = out.PatchOrder';
            tmp.env = repelem(blockI,numel(tmp.leaveT))';

            simData = [simData; tmp]; % concatenate onto previous block
        end

        simSubjData = PrepSubjData(simData, block, Env, blockFlag); % transform into state actions ready for fitting

        % fit all models back to the simulated data
        [BIC, iBEST, BEST] = fitAll(allModels, simSubjData, Env);
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

xticklabels(modelNames)
yticklabels(modelNames)
end


