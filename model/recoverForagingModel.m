%% parameter recovery

clear
close all;
addpath(genpath('./foraging/code'))

%% user options
model = 5; % model type - see model table to check number to choose
blockFlag = 'combined'; %% either 'combined' (fit as one continuous task) or 'separate' (fit rich and poor as separate blocks)
% note that if choosing 'separate', select below which block to look at
% (rich 1 or poor 2) 
block = 1;

nSim = 250; % number of iterations to recover

%% set up agent and task
SetUpEnviron % get the environment parameters
blockOrder = [1 2]; 

[~, ~, paramsIndex] = buildForagingModel(model);
numParams = sum(paramsIndex);

blockNames = {'rich', 'poor'};
paramNames = {'alpha patch/Q', 'alpha rho', 'beta'};
paramNames = paramNames(paramsIndex);

%% fmincon options
lowerBounds = [0, 0, 0];lowerBounds = lowerBounds(paramsIndex);
upperBounds = [1, 1, 100];upperBounds = upperBounds(paramsIndex);
options = optimoptions('fmincon','Display','none'); % don't display

nStarts = 5; % how many starting locations to begin fitting (avoid local minima)

minNLLFitParams = zeros([nSim, numParams]);
simParams = zeros([nSim, numParams]);

for ii = 1:nSim
    ii

    %% simulate data
    simData = table; % container table for simulated data

    % agent parameters - same parameters for both blocks
    params = [defineBoundedParam(0, 1), defineBoundedParam(0, 0.05), exprnd(1)];
    params = params(paramsIndex);
    simParams(ii,:) = params;

    for blockI = 1:2

        blockType = blockOrder(blockI); % for recovery, assume always seen rich block first (same order as blockI)
        Env.BlockPatchOrder = Env.PatchOrder{blockType}; % set the current patch order

        if blockI == 2 && contains(blockFlag, 'combined') % we need to specify start values if it's the 2nd block and we're doing combined fitting of both blocks together
            startValues.Rho = out.Rho(end); % first value is last value from previous block
            startValues.EstimatedPatchRR = out.EstimatedPatchRR(end);
        else
            startValues.Rho = 0;
            startValues.EstimatedPatchRR = 0;
        end

        sim_f = buildForagingModel(model, Env, [], startValues); % get function handle for this model
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

    %% fit back to the simulated data 

    [~, NLL_f] = buildForagingModel(model, Env, simSubjData); % need to update NLL function with new subject data each time

    NLLEval = zeros([nStarts, 1]);
    fitParams = zeros([nStarts, numParams]);

    % recover data
    parfor mm = 1:nStarts
        params0 = [defineBoundedParam(0, 1), defineBoundedParam(0, 1), exprnd(1)];
        params0 = params0(paramsIndex);
        [fitParams(mm,:),NLLEval(mm)] = fmincon(NLL_f,params0,[],[],[],[],lowerBounds, upperBounds,[],options);
    end

    minNLL = min(NLLEval);   % minimum negative log likelihood over all starting positions
    ix = find(minNLL == NLLEval);    % indices of location of minimum, to find the corresponding best fit parameters
    minNLLFitParams(ii,:) = fitParams(ix(1),:); % get corresponding parameter values at lowest NLL - take 1st index if two lowest NLL

end

%% plot correlation - simulated value versus actual value
close all
figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');

% find 'bad' alpha rho values
%thresh = 0.25;
%ind = abs(simParams(:,2) - minNLLFitParams(:,2)) > thresh;

% for each parameter, plot sim vs fit
for i= 1:numParams
    ax = nexttile;
    plot(simParams(:,i), minNLLFitParams(:,i), 'o', 'markersize', 8, 'linewidth', 1); hold on;
    %plot(simParams(ind, i), minNLLFitParams(ind,i), 'o', 'color', 'blue', 'markersize', 8, 'linewidth', 1,'markerfacecolor', [1 1 1]*0.5)
    title(paramNames{i})
    xlabel(sprintf('Simulated %s', paramNames{i}))
    ylabel(sprintf('Fit %s', paramNames{i}))
    if contains(paramNames{i},'beta')
        corrcoef(log(simParams(:,i)),log(minNLLFitParams(:,i))) % do correlation between logs of beta's
        set(ax,'yscale' ,'log'); set(ax,'xscale' ,'log')
    else
        corrcoef(simParams(:,i),minNLLFitParams(:,i))
    end
end
title(tl, 'Parameter recovery')

% plot trade off between parameters

figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');
ax = nexttile;
plot(minNLLFitParams(:,1), minNLLFitParams(:,2), 'o', 'markersize', 8, 'linewidth', 1)
title('fit alpha patch/Q vs alpha rho')
xlabel('fit alpha patch/Q')
ylabel('fit alpha rho')
set(ax,'yscale' ,'log')

ax = nexttile;
plot(minNLLFitParams(:,1), minNLLFitParams(:,3), 'o', 'markersize', 8, 'linewidth', 1)
title('fit alpha patch/Q vs fit beta')
xlabel('fit alpha patch/Q')
ylabel('fit beta')

ax = nexttile;
plot(minNLLFitParams(:,2), minNLLFitParams(:,3), 'o', 'markersize', 8, 'linewidth', 1)
title('fit alpha rho vs fit beta')
xlabel('fit alpha rho')
ylabel('fit beta')
set(ax,'yscale' ,'log')

%% checks - simulation and fitting scripts end up with the same outputs

%  TestParams = [0.4, 0.1, 1];
% % 
% % % simulate data
% block = 1; 
%  out = simulate_MVT_full_RW_softmax(TestParams, block, Env); % change this depending on model to test
%  [SubjNLLEval, RecoveredData] = NLL_MVT_full_RW_softmax(TestParams, Env, out);
% out.PAction == RecoveredData.PAction
% out.Q == RecoveredData.Q
% out.Rho == RecoveredData.Rho
% out.NextPatch == RecoveredData.NextPatch
