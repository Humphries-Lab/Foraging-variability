function [agent] = loadAgents(study,modelOptions,params,nSim,task,fitDataFlag)

rich = 1; poor = 2;

% select parameters depending on model

if fitDataFlag == 0 % if simulating from scratch

    agent.summary.blockOrder = [repmat([1 2], [nSim/2, 1]); repmat([2 1], [nSim/2, 1])];
    agent.summary.experiencedAvgRR = zeros([nSim,2]);

    agent.params = cell([1 2]);
    agent.data = cell([1 2]);

    if strcmp(modelOptions.blockPresentation,'combined')
        params.poor = params.rich; % only read in one set of parameters
    end

    T = array2table(repmat(params.rich,[nSim,1])); T.Properties.VariableNames = params.names; tmp{rich} = T;
    T = array2table(repmat(params.poor,[nSim,1])); T.Properties.VariableNames = params.names; tmp{poor} = T; clear T

    for iB = 1:numel(task.environNames) 

        % set up action policy parameters
        switch modelOptions.actionPolicy

            case 'softmax'

                switch modelOptions.betaFunction
                    case 'fit'
                        agent.params{iB}.beta = tmp{iB}.beta;
                    case 'linear' | 'exponential'
                        agent.params{iB}.lambda = tmp{iB}.lambda;
                    case 'twoLinear' | 'twoExponential'
                        agent.params{iB}.lambda = tmp{iB}.lambda;
                        agent.params{iB}.gamma = tmp{iB}.gamma;
                end

            case 'epsilon-greedy'
                agent.params{iB}.epsilon = tmp{iB}.epsilon;

            case 'epsilon-softmax'
                agent.params{iB}.beta = tmp{iB}.beta;
                agent.params{iB}.epsilon = tmp{iB}.epsilon;
        end

        % set up learning rates for average RR
        if modelOptions.learnRho
            agent.params{iB}.alphaRho = tmp{iB}.alphaRho;
        end

        % set up learning rates for patch RR
        if modelOptions.learnPatchRR
            agent.params{iB}.alphaPatch = tmp{iB}.alphaPatch;
        end

        agent.params{iB} = struct2table(agent.params{iB});

        % set up empty container to store their data
        agent.data{iB}.pAction = zeros(task.blockTime,2); % probability of selecting leave or stay
        agent.data{iB}.action = zeros(task.blockTime,1); % what action taken
        agent.data{iB}.reward = zeros(task.blockTime,1); % the reward earned in each state in the block
        agent.data{iB}.estPatchRR = zeros(task.blockTime,1);  % estimated patchRR
        agent.data{iB}.rho = zeros(task.blockTime,1); % estimated averageRR
        
        agent.data{iB} = struct2table(agent.data{iB});

    end


else

    switch study
        case 'leheron'
            load(sprintf('~/Dropbox/foraging/outputs/fitting/fitting_results_%s_M%d', modelOptions.blockPresentation, modelOptions.modelNumber), 'minNLLFitParams');
            load('~/Dropbox/foraging/raw_data/BlockOrder.mat'); % load order of block presentation
            load('~/Dropbox/foraging/raw_data/experiencedAvgRR.mat'); % load their real average RR in the 2 blocks
            agent.blockOrder = BlockOrder;
            agent.experiencedAvgRR = experiencedAvgRR;
            agent.ID = [1:size(minNLLFitParams,1)]';
            %             if contains(blockPresentation, 'separate')
            %                 sim = minNLLFitParams;
            %             elseif contains(blockPresentation, 'combined')
            %                 sim = cat(3, minNLLFitParams, minNLLFitParams); % duplicate for aligned block indexing later
            %             end
    end

end

agent.summary = struct2table(agent.summary);
