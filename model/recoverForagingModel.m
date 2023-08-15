%% ---- Script to recover MVT model for patch foraging tasks ---- %%
% Emma Scholey
% latest update 14 August 2023

clear
close all

addpath('./helperFunctions')
addpath('./analyticalComputation')

%% user options

% study options
study = 'leheron'; % study to simulate/fit data to.

% model options
modelNum = 1; % model type - see model table to check number to choose
modelTable = readtable('./foragingModelTable.xlsx'); 

% simulation options
funcOptions.type = 'simulate_new'; % 'simulate_new' if simulating new parameters, 'simulate_fit' if simulating already fit parameters for each subject
funcOptions.blockPresentation = 'separate'; % either 'combined' (fit as one continuous task) or 'separate' (fit rich and poor as separate blocks)
funcOptions.nSim = 100;

%% set up

% load task
task = buildTask(study,funcOptions.blockPresentation); % set up task

% load dataframe container for simulations
allData = buildData(study,task,funcOptions);

% load model
model = table2struct(modelTable(modelTable.modelNumber == modelNum,:));

% load agent parameters 
funcOptions.type = 'fit'; % we want to generate new parameters from scratch, not load existing. Temporarily change flag.
allParams = buildParams(study,task,model,funcOptions); clear params
model.paramNames = allParams.names;
funcOptions.type = 'simulate_new'; % revert back 

%% run simulations

for iS = 1:funcOptions.nSim

    for iB = 1:task.nEnviron

        iS

        % what is the agent data for this block?
        agent.currentBlock = allData.blockOrder(iS,iB); % get current block for this subject based on their block order

        agent.experiencedAvgRR = allData.experiencedAvgRR(iS,:); 
        agent.nStates = allData.nStates(iS,agent.currentBlock);
        agent.data = allData.data{agent.currentBlock}; % empty generic container
        agent.patchOrder = task.patchOrder{agent.currentBlock}'; % generic patch order for all participants

        agentParams = table2array(allParams.params(iS,:)); % specific for the subject. Use same parameter combo for each block

        if iB == 2 && strcmp(funcOptions.blockPresentation, 'combined') % specify start values if it's the 2nd block and we're doing combined fitting of both blocks together
            agent.data.rho(1) = out.rho(end);
        end

        [~,results,out] = simulate_MVT_model(task,model,agent,agentParams,funcOptions);

        simData.data{iS}{iB} = out; 
        simData.results{iS}{iB} = results;
        simData.nObservations(iS,iB) = sum(out.action == 2);
    end
end

%% fit back to the simulated data 
funcOptions.type = 'fit'; % 'simulate_new' if simulating new parameters, 'simulate_fit' if simulating already fit parameters for each subject
funcOptions.nSim = 5; % how many search points 

% load random set of start parameters for fmincon search
searchParams = buildParams(study,task,model,funcOptions);
paramArray = table2array(allParams.params);

options = optimoptions('fmincon','Display','none'); % don't display
lb = allParams.lb;
ub = allParams.ub;

% initialise containers
minNLL = zeros([numel(simData.data), task.blockNum]);
minNLLFitParams = zeros([numel(simData.data), allParams.nParams, task.blockNum]);
BIC = zeros([numel(simData.data), task.blockNum]);
AIC = zeros([numel(simData.data), task.blockNum]);

for iS = 1:numel(simData.data) % for all simulated combinations

    subj.experiencedAvgRR = allData.experiencedAvgRR(iS,:); 

    for iB = 1:task.blockNum
        iS
        subj.nStates = allData.nStates(iS,iB);
        subj.data = simData.data{iS}{iB};
        subj.data.estPatchRR = zeros(subj.nStates+1,1);   % refresh for recovery
        subj.data.rho = zeros(subj.nStates+1,1);    % refresh for recovery
        subj.patchOrder = simData.results{iS}{iB}.patchOrder;
        subj.env = simData.results{iS}{iB}.env;

        NLLEval = zeros([funcOptions.nSim, 1]);
        FitParams = zeros([funcOptions.nSim, allParams.nParams]);

        % Run fmincon
        for ii = 1:funcOptions.nSim
            params0 =  paramArray(ii,:);

            f = @(x0)simulate_MVT_model(task,model,subj,x0,funcOptions);
            [FitParams(ii,:),NLLEval(ii)] = fmincon(f,params0,[],[],[],[],lb,ub,[],options);
        end

        % Find the best fitting parameter values
        minNLL = min(NLLEval);   % minimum negative log likelihood over all starting positions
        ix = find(minNLL == NLLEval);    % indices of location of minimum, to find the corresponding best fit parameters
        minNLLFitParams(iS,:,iB) = FitParams(ix(1),:); % get corresponding parameter values at lowest NLL

        % Calculate BIC/AIC
        BIC(iS,iB) = allParams.nParams * log(simData.nObservations(iS,iB)) + 2*minNLL;
        AIC(iS,iB) = 2/simData.nObservations(iS,iB) * minNLL + 2 * allParams.nParams/simData.nObservations(iS,iB);

    end
end
clear FitParams NLLEval iS ix


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
