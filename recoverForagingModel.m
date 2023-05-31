%% parameter recovery

clear
close all;
addpath(genpath('~/Dropbox/foraging/code'))

model = 25;

block = 2; % doesn't really matter
SetUpEnviron % get the environment parameters
blockNames = {'rich', 'poor'};
paramNames = {'alpha patch/Q', 'alpha rho', 'beta'};

%% simulation and recovery of parameters

lowerBounds = [0, 0, 0];
upperBounds = [1, 1, 100];
options = optimoptions('fmincon','Display','none'); % don't display

nSim = 250; % number of iterations to recover
nStarts = 5; % how many starting locations to begin fitting (avoid local minima)

numParams = length(lowerBounds);

minNLLFitParams = zeros([nSim, numParams]);
simParams = zeros([nSim, numParams]);

for ii = 1:nSim
    ii
    % agent parameters
    simParams(ii,:) = [defineBoundedParam(0, 1), defineBoundedParam(0, 1), exprnd(10)];

    % simulate data

    %out = simulateM1_MVT(simParams(ii,:), block, Env); % change this depending on model to test
    %out = simulateM2_MVT_RW(simParams(ii,:), block, Env); % change this depending on model to test
    %out = simulateM3_RLOn(simParams(ii,:), block, Env); % change this depending on model to test
    %out = simulateM4_RLOn(simParams(ii,:), block, Env); % change this depending on model to test
    %out = simulateM5_RLOn(simParams(ii,:), block, Env); % change this depending on model to test
    %out = simulateM21_RLOff(simParams(ii,:), block, Env); % change this depending on model to test
    out = simulateM25_RLOn(simParams(ii,:), block, Env); % change this depending on model to test
    %out = simulateM26_RLOn(simParams(ii,:), block, Env); % change this depending on model to test

    out.PatchOrder = Env.PatchOrder{block}; % give unlimited patch order so that fitting back doesn't break 

    NLLEval = zeros([nStarts, 1]);
    fitParams = zeros([nStarts, numParams]);

    % recover data
    parfor mm = 1:nStarts
        params0 = [defineBoundedParam(0, 1), defineBoundedParam(0, 1), exprnd(10)];

        %f = @(x)NLL_M1_MVT(x, Env, out, block);  % need to tell MVT proper which block, to get optimal average RR as threshold
        %f = @(x)NLL_M2_MVT_RW(x, Env, out);
        %f = @(x)NLL_M3_RLOn(x, Env, out);
        %f = @(x)NLL_M4_RLOn(x, Env, out);
        %f = @(x)NLL_M5_RLOn(x, Env, out);

        %f = @(x)NLL_M21_RLOff(x, Env, out);
        f = @(x)NLL_M25_RLOn(x, Env, out);
        %f = @(x)NLL_M26_RLOn(x, Env, out);

        [fitParams(mm,:),NLLEval(mm)] = fmincon(f,params0,[],[],[],[],lowerBounds, upperBounds,[],options);
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

ax = nexttile;
plot(minNLLFitParams(:,1), minNLLFitParams(:,3), 'o', 'markersize', 8, 'linewidth', 1)
title('fit alpha patch/Q vs fit beta')
xlabel('fit alpha patch/Q')
ylabel('fit beta')
set(ax,'yscale' ,'log')

ax = nexttile;
plot(minNLLFitParams(:,2), minNLLFitParams(:,3), 'o', 'markersize', 8, 'linewidth', 1)
title('fit alpha rho vs fit beta')
xlabel('fit alpha rho')
ylabel('fit beta')
set(ax,'yscale' ,'log')

%% checks - simulation and fitting scripts end up with the same outputs

TestParams = [0.4, 0.1, 1];

% simulate data
out = simulateM21_RLOff(TestParams, block, Env); % change this depending on model to test
[SubjNLLEval, RecoveredData] = NLL_M21_RLOff(TestParams, Env, out);
out.PAction == RecoveredData.PAction
out.Q == RecoveredData.Q
out.Rho == RecoveredData.Rho
out.NextPatch == RecoveredData.NextPatch
