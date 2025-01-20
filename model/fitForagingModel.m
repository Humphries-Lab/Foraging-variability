%% fitting models to foraging data
% Emma Scholey 9 Jun 2022
% latest update 4 August 2023

clear
close all

addpath('./helperFunctions')

%% user options

% study options
fitOptions.study = 'contrerashuerta'; % Options are leheron (field-human), contrerashuerta (berry-human), or kane (rat)

% model options
modelTable = readtable('./foragingModelTable.xlsx');

% fitting options
fitOptions.type = 'fit';
fitOptions.nSim = 50; % how many starts/iterations for fmincon search

%% set up model and task

% load task
task = buildTask(fitOptions.study);

% load participant data
allData = buildData(task,fitOptions);
nSub = size(allData.data,2);

for modelNum = [17] % for all models
    % load model
    model = table2struct(modelTable(modelTable.modelNumber == modelNum,:));

    % load random set of start parameters for fmincon search
    allParams = buildParams(model,fitOptions);
    model.paramNames = allParams.names;

    %% fitting for each person in group with different starting points
    options = optimoptions('fmincon','Display','none'); % don't display
    lb = allParams.lb; % define lower and upper bounds, lb, ub
    ub = allParams.ub;

    % initialise containers
    minNLL = zeros([nSub 1]);
    minNLLFitParams = zeros([nSub, allParams.nParams]);
    BIC = zeros([nSub, 1]);
    AIC = zeros([nSub, 1]);

    paramArray = table2array(allParams.params);

    for iS = 1:nSub

        subj.experiencedAvgRR = allData.experiencedAvgRR(iS,:); % real experienced avgRR

        iS

        subj.nStates = allData.nStates(iS);
        subj.data = allData.data{iS};
        if model.learnRho == 1
            subj.data.rho(1) = mean(task.optAvgRR); % initialise rho to the mean optimal avg reward rate, for learning models only
        end
        subj.patchOrder = allData.patchOrder{iS};
        subj.env = allData.env{iS};
        subj.switchIndex = allData.switchIndex{iS};

        NLLEval = zeros([fitOptions.nSim, 1]);
        FitParams = zeros([fitOptions.nSim, allParams.nParams]);

        % Run fmincon
        parfor ii = 1:fitOptions.nSim
            params0 =  paramArray(ii,:);

            % test code (check negLL computed correctly)
            %[NLL, results, out] = simulate_MVT_model(task,model,subj,[0.31, 0.0076,-1.99],fitOptions);

            f = @(x0)simulate_MVT_model(task,model,subj,x0,fitOptions);
            [FitParams(ii,:),NLLEval(ii)] = fmincon(f,params0,[],[],[],[],lb,ub,[],options);
        end

        % Find the best fitting parameter values
        minNLL = min(NLLEval);   % minimum negative log likelihood over all starting positions
        ix = find(minNLL == NLLEval);    % indices of location of minimum, to find the corresponding best fit parameters
        minNLLFitParams(iS,:) = FitParams(ix(1),:); % get corresponding parameter values at lowest NLL

        % Calculate BIC/AIC
        BIC(iS) = allParams.nParams * log(allData.nObservations(iS)) + 2*minNLL;
        AIC(iS) = 2/allData.nObservations(iS) * minNLL + 2 * allParams.nParams/allData.nObservations(iS);

    end
    clear FitParams NLLEval iS ix

    %% save results
    m = median(minNLLFitParams);

    minNLLFitParams = array2table(minNLLFitParams, "VariableNames",model.paramNames);
    % save_name = sprintf('../data/fitting_data/fitting_results_M%d_%s', model.modelNumber,fitOptions.study);
    % save(save_name, 'AIC', 'BIC', 'minNLLFitParams')
end
