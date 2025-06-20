%% ---- Script to recover foraging model for patch foraging task - model identifiability ---- %%
% Emma Scholey
% created 17 April 2025
% latest update 17 April 2025

clear
close all
addpath '/rds/projects/a/appsmaj-motivation-social-neuro/Emma/spm'
addpath('./helperFunctions')

%% user options

sim_iter = 100; 
fit_iter = 5;

% study options
funcOptions.study = 'contrerashuerta'; % Options are leheron (field-human), contrerashuerta (berry-human), or kane (rat)

% model options
modelNum = [1:5]; % model type - see model table to check number to choose

modelTable = readtable('./foragingModelTable.xlsx');

recoveredParams = nan(sim_iter, 4, length(modelNum), length(modelNum)); % 4 = max nParams

%% set up

for iM = 2%1:numel(modelNum)

    % simulation options
    funcOptions.type = 'simulate_new'; % 'simulate_new' if simulating new parameters, 'simulate_fit' if simulating already fit parameters for each subject
    funcOptions.nSim = sim_iter;
    % load task
    task = buildTask(funcOptions.study); % set up task

    % load dataframe container for simulations
    allData = buildData(task, funcOptions);

    % load model
    model = table2struct(modelTable(modelTable.modelNumber == modelNum(iM),:));

    % load agent parameters
    funcOptions.type = 'recover'; % we want to generate new parameters from scratch, not load existing. Temporarily change flag.
    allParams = buildParams(model,funcOptions); clear params
    model.paramNames = allParams.names;

    simulatedParams = allParams.params;
    funcOptions.type = 'simulate_new'; % revert back

    %% run simulations

    simData = cell(1,funcOptions.nSim);
    simLT = cell(1,funcOptions.nSim);

    for iS = 1:funcOptions.nSim

        agent.experiencedAvgRR = allData.experiencedAvgRR(iS,:); % real experienced avgRR (calculated post-hoc)
        initialiseRho = mean(task.optAvgRR); % set up rho for block 1 (only needed for learning models)

        estPatchRR = []; 
        rho = [];
        action = [];
        patchOrder = [];
        env = [];
        switchIndex = [];

        for iB = 1:task.nBlocks % why are we doing by block? because we don't know how many trials (patches) in each time-restricted block, so we need to simulate each block separately
            agent.nStates = allData.nStates;
            agent.data = allData.data; % empty generic container
            agent.data.rho(1) = initialiseRho; % carry over rho from previous block

            agent.currentEnv = allData.blockOrder(iS,iB);

            agent.patchOrder = task.patchOrder{agent.currentEnv}'; % generic patch order for all participants

            agentParams = allParams.params{iS,:}; % specific for the subject, depending on environment

            [~,lt,results] = simulate_MVT_model(task,model,agent,agentParams,funcOptions); % run the simulation

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
        nObservations(iS) = sum(action == 2);
    end

    %% fit all models back to the simulated data

    for jM = 1:numel(modelNum)

        model = table2struct(modelTable(modelTable.modelNumber == modelNum(jM),:));
        funcOptions.type = 'fit'; % 'simulate_new' if simulating new parameters, 'simulate_fit' if simulating already fit parameters for each subject
        funcOptions.nSim = fit_iter; % how many search points

        % load random set of start parameters for fmincon search
        searchParams = buildParams(model,funcOptions);
        model.paramNames = searchParams.names;
        paramArray = table2array(searchParams.params);

        options = optimoptions('fmincon','Display','none'); % don't display
        lb = searchParams.lb;
        ub = searchParams.ub;

        for iS = 1:allData.nSubj

            disp(num2str([modelNum(iM),modelNum(jM), iS]))

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
            BIC(iS,jM,iM) = searchParams.nParams * log(nObservations(iS)) + 2*minNLL;
            ix = find(minNLL == NLLEval);    % indices of location of minimum, to find the corresponding best fit parameters
            bestP = FitParams(ix(1),:);
            Pj = numel(bestP);
            recoveredParams(iS, 1:Pj, jM, iM) = bestP;

            minNLL_all(iS,jM,iM) = minNLL;
        end
    end 

    posteriorProbabilities = BICposterior(BIC(:,:,iM));
    [~, ~, exceedance_probabilities(iM,:), protected_XP(iM,:)] = spm_BMS(posteriorProbabilities);

end
modelNames = {'vary \beta', 'vary \beta vary c', 'vary \beta fix c', 'fix \beta vary c', 'fix \beta fix c'};

save(sprintf('../data/simulation_data/XP_model_identifiability_%s_v4', funcOptions.study), 'exceedance_probabilities','protected_XP', 'modelNames', 'allParams', 'recoveredParams', 'posteriorProbabilities')

%% plot heatmap

figure
h = heatmap(exceedance_probabilities,'MissingDataColor','w', 'ColorMap', sky, 'GridVisible', 'off','CellLabelColor','none'); % remember to transpose 
labels = modelNames;
h.XDisplayLabels = labels;
h.YDisplayLabels = labels;
h.XLabel = 'estimated';
h.YLabel = 'simulated';
set(gca,'FontSize',16)
