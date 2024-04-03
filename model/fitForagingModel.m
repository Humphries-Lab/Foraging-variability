%% fitting models to foraging data
% Emma Scholey 9 Jun 2022
% latest update 4 August 2023

clear
close all

addpath('./helperFunctions')
addpath('./analyticalComputation')

%% user options

% study options
fitOptions.study = 'leheron'; % study to simulate/fit data to.

% model options
modelTable = readtable('./foragingModelTable.xlsx');

% fitting options
fitOptions.type = 'fit'; % not simulating data here
fitOptions.nSim = 8; % how many starts/iterations for fmincon search

%% set up model and task

% load task
task = buildTask(fitOptions.study);

% load participant data
allData = buildData(task,fitOptions);
nSub = size(allData.data,2);

for modelNum = 6:8 % for all models
    % load model
    model = table2struct(modelTable(modelTable.modelNumber == modelNum,:));

    % load random set of start parameters for fmincon search
    allParams = buildParams(task,model,fitOptions);
    model.paramNames = allParams.names;

    %% fitting for each person in group with different starting points
    options = optimoptions('fmincon','Display','none'); % don't display
    lb = allParams.lb;
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


    %% plots

    %     close all
    %
    %     if strcmp(fitOptions.blockPresentation,'separate')
    %         plotNames = task.environNames;
    %     elseif strcmp(fitOptions.blockPresentation,'combined')
    %         plotNames = {'combined'};
    %     end
    %
    %     if allParams.nParams > 1
    %         combinations = nchoosek(1:allParams.nParams,2);
    %
    %         for iB = 1:task.numFitBlocks
    %             figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');
    %
    %             for i= 1:size(combinations,1)
    %                 nexttile;
    %                 scatter(minNLLFitParams(:,combinations(i,1),iB),minNLLFitParams(:,combinations(i,2),iB))
    %                 xlabel(sprintf('Fit %s', model.paramNames{:,combinations(i,1)}))
    %                 ylabel(sprintf('Fit %s', model.paramNames{:,combinations(i,2)}))
    %                 title(plotNames{iB})
    %             end
    %         end
    %         title(tl, 'Best fit parameter distributions')
    %     else
    %         figure
    %         tl = tiledlayout('flow', 'TileSpacing', 'Compact');
    %         for iB = 1:task.numFitBlocks
    %             nexttile;
    %             boxchart(minNLLFitParams(:,:,iB))
    %             xlabel(sprintf('Fit %s', model.paramNames{1}))
    %             title(plotNames{iB})
    %         end
    %     end


    %% save results
    m = median(minNLLFitParams);

    minNLLFitParams = array2table(minNLLFitParams, "VariableNames",model.paramNames);
    save_name = sprintf('../data/fitting_data_240117/fitting_results_M%d_%s', model.modelNumber,fitOptions.study);
    save(save_name, 'AIC', 'BIC', 'minNLLFitParams')
end
