%% ---- Script to recover AET model for predator-prey task ---- %%
% Emma Scholey
% latest update 14 August 2023

clear
close all

addpath('./helperFunctions')

%% user options

funcOptions.study = 'kane';

% model options
modelNum = 3; % model type - see model table to check number to choose
modelTable = readtable('./foragingModelTable.xlsx'); 


%% set up

% simulation options
funcOptions.type = 'simulate_new'; % 'simulate_new' if simulating new parameters, 'simulate_fit' if simulating already fit parameters for each subject
funcOptions.nSim = 100;

% load task
task = buildTask(funcOptions.study); % set up task

% load dataframe container for simulations
allData = buildData(task, funcOptions);

% load model
model = table2struct(modelTable(modelTable.modelNumber == modelNum,:));

% load agent parameters 
funcOptions.type = 'recover'; % we want to generate new parameters from scratch, not load existing. Temporarily change flag.
allParams = buildParams(model,funcOptions); clear params
model.paramNames = allParams.names;

funcOptions.type = 'simulate_new'; % revert back 

%% run simulations

for iS = 1:funcOptions.nSim

    iS 
    agent.experiencedAvgRR = allData.experiencedAvgRR(iS,:); % real experienced avgRR (calculated post-hoc)
    initialiseRho = mean(task.optAvgRR); % set up rho for block 1 (only needed for learning models)

    estPatchRR = [];
    rho = [];
    action = [];
    patchOrder = [];
    env = [];

    for iB = 1:task.nBlocks % why are we doing by block? because we don't know how many trials (patches) in each time-restricted block, so we need to simulate each block separately
        
        agent.nStates = allData.nStates;
        agent.data = allData.data; % empty generic container
        agent.data.rho(1) = initialiseRho; % carry over rho from previous block

        agent.currentEnv = allData.blockOrder(iS,iB);

        agent.patchOrder = task.patchOrder{agent.currentEnv}'; % generic patch order for all participants

        agentParams = allParams.params{iS,:}; % specific for the subject, depending on environment

        [NLL,lt,results] = simulate_MVT_model(task,model,agent,agentParams,funcOptions); % run the simulation

        initialiseRho = results.rho(end); % what was the last rho from previous block (to carry over)

        % Concatenate blocks together
        estPatchRR = [estPatchRR; results.estPatchRR];
        rho = [rho; results.rho];
        action = [action; results.action];

        patchOrder = [patchOrder; lt.patchOrder];
        env = [env; lt.env];


    end

    patchOrder = [patchOrder; 2]; % add a last patch to account for fact that last patch is left early, so not logged in simulated data
    env = [env; env(end)];

    switchIndex = diff(env) ~= 0; switchIndex = [1; switchIndex];
    simLT{iS} = struct('patchOrder', patchOrder, 'env', env, 'switchIndex', switchIndex);
    simData{iS} = struct('estPatchRR', estPatchRR, 'rho', rho, 'action', action);
end

%% fit back to the simulated data 
funcOptions.type = 'fit'; % 'simulate_new' if simulating new parameters, 'simulate_fit' if simulating already fit parameters for each subject
funcOptions.nSim = 5; % how many search points 

% load random set of start parameters for fmincon search
searchParams = buildParams(model,funcOptions);
paramArray = table2array(searchParams.params);

options = optimoptions('fmincon','Display','none'); % don't display
lb = searchParams.lb;
ub = searchParams.ub;

% initialise containers
minNLLFitParams_recovered = zeros([allData.nSubj allParams.nParams]);

for iS = 1:allData.nSubj

    iS

    % Extract each subject
    subj.data = simData{iS};
    subj.experiencedAvgRR = allData.experiencedAvgRR(iS,:); % real experienced avgRR
    
    subj.nStates = allData.nStates*task.nBlocks;
    if model.learnRho == 1
        subj.data.rho(1) = mean(task.optAvgRR); % initialise rho to the mean optimal avg reward rate, for learning models only
    end
    subj.patchOrder = simLT{iS}.patchOrder;
    subj.env = simLT{iS}.env;
    subj.switchIndex = simLT{iS}.switchIndex;

    NLLEval = zeros([funcOptions.nSim, 1]);
    FitParams = zeros([funcOptions.nSim, searchParams.nParams]);

    % Run fmincon
    parfor ii = 1:funcOptions.nSim
        params0 =  paramArray(ii,:);

        f = @(x0)simulate_MVT_model(task,model,subj,x0, funcOptions);
        [FitParams(ii,:),NLLEval(ii)] = fmincon(f,params0,[],[],[],[],lb,ub,[],options);
    end

    % Find the best fitting parameter values
    minNLL = min(NLLEval);   % minimum negative log likelihood over all starting positions
    ix = find(minNLL == NLLEval);    % indices of location of minimum, to find the corresponding best fit parameters
    minNLLFitParams_recovered(iS,:) = FitParams(ix(1),:); % get corresponding parameter values at lowest NLL

end

%% plots
close all

figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');

for i= 1:allParams.nParams
    nexttile;
    scatter(allParams.params{:,i},minNLLFitParams_recovered(:,i))
    
    xlabel(['Simulated ' , model.paramNames{i}])
    ylabel(['Fit ' , model.paramNames{i}])
    model.paramNames{i}

    corr([allParams.params{:,i}, minNLLFitParams_recovered(:,i)], 'type', 'Spearman')
end


% plot trade off between parameters 
if allParams.nParams > 1
    combinations = nchoosek(1:allParams.nParams,2);
    figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');

    for i= 1:size(combinations,1)
        nexttile;
        scatter(minNLLFitParams_recovered(:,combinations(i,1)),minNLLFitParams_recovered(:,combinations(i,2)))
        xlabel(sprintf('Fit %s', model.paramNames{:,combinations(i,1)}))
        ylabel(sprintf('Fit %s', model.paramNames{:,combinations(i,2)}))

       disp(model.paramNames{:,combinations(i,1)})
       disp(model.paramNames{:,combinations(i,2)})
       [r, p] = corr(minNLLFitParams_recovered(:,combinations(i,1)),minNLLFitParams_recovered(:,combinations(i,2)), 'type', 'Spearman')
    end
end

params.simulated = table2array(allParams.params);
params.estimated = minNLLFitParams_recovered; 
params.names = {'\beta rich', '\beta poor', 'c'};

save(sprintf('../data/simulation_data/M%d_parameter_recovery_%s', modelNum, funcOptions.study), 'params')

%% plot heatmap

r = corr(table2array(allParams.params),minNLLFitParams_recovered, type="Spearman");
r = round(r,2);
figure
h = heatmap(r,'MissingDataColor','w', 'ColorMap', sky, 'GridVisible', 'off');
labels = ["\beta rich","\beta poor","c"];
h.XDisplayLabels = labels;
h.YDisplayLabels = labels; 
h.XLabel = 'simulated';
h.YLabel = 'estimated';

%print([sprintf('../../figures/panels/parameter_recovery_M%d_',modelNum), funcOptions.version],'-dsvg')
