function [agentParams] = buildParams(study,task,modelOptions,funcOptions)

switch funcOptions.type
    case 'simulate_new'

        if strcmp(funcOptions.blockPresentation,'combined')
            funcOptions.params.poor = funcOptions.params.rich; % only read in one set of parameters
        end

        % convert rich parameters to table to index by name 
        T = array2table(repmat(funcOptions.params.rich,[funcOptions.nSim,1])); 
        T.Properties.VariableNames = {'beta', 'epsilon', 'alphaRho', 'alphaPatch', 'lambda', 'gamma'}; tmp{1} = T;
        
        % convert poor parameters to table to index by name 
        T = array2table(repmat(funcOptions.params.poor,[funcOptions.nSim,1])); 
        T.Properties.VariableNames = {'beta', 'epsilon', 'alphaRho', 'alphaPatch', 'lambda', 'gamma'}; tmp{2} = T; clear T

        for iB = 1:numel(task.environNames)

            % set up action policy parameters
            switch modelOptions.actionPolicy

                case 'softmax'

                    switch modelOptions.betaFunction
                        case 'fit'
                            agentParams.params{iB}.beta = tmp{iB}.beta;
                        case {'linear', 'exponential'}
                            agentParams.params{iB}.lambda = tmp{iB}.lambda;
                        case {'twoLinear', 'twoExponential'}
                            agentParams.params{iB}.lambda = tmp{iB}.lambda;
                            agentParams.params{iB}.gamma = tmp{iB}.gamma;
                    end

                case 'epsilon-greedy'
                    agentParams.params{iB}.epsilon = tmp{iB}.epsilon;

                case 'epsilon-softmax'
                    agentParams.params{iB}.beta = tmp{iB}.beta;
                    agentParams.params{iB}.epsilon = tmp{iB}.epsilon;
            end

            % set up learning rates for average RR
            if modelOptions.learnRho || strcmp(modelOptions.rhoFunction, 'rw')
                agentParams.params{iB}.alphaRho = tmp{iB}.alphaRho;
            end

            % set up learning rates for patch RR
            if modelOptions.learnPatchRR
                agentParams.params{iB}.alphaPatch = tmp{iB}.alphaPatch;
            end

            agentParams.params{iB} = struct2table(agentParams.params{iB});

        end
            agentParams.names = agentParams.params{1}.Properties.VariableNames; 

    case 'simulate_fit'

        switch study
            case 'leheron'
                load(sprintf('~/Dropbox/foraging/outputs/fitting/fitting_results_%s_M%d', funcOptions.blockPresentation, modelOptions.modelNumber), 'minNLLFitParams');
                if contains(funcOptions.blockPresentation, 'separate')
                    agentParams = minNLLFitParams;
                elseif contains(funcOptions.blockPresentation, 'combined')
                    agentParams = cat(3, minNLLFitParams, minNLLFitParams); % duplicate for aligned block indexing later
                end
            case 'contrerashuerta'
            case 'kane'
        end

    case 'fit'
        % generate set of parameters to start fmincon search

        % set up action policy parameters
        switch modelOptions.actionPolicy

            case 'softmax'

                switch modelOptions.betaFunction
                    case 'fit'
                        agentParams.params.beta = exprnd(3,[funcOptions.nStarts,1]);
                    case {'linear', 'exponential'}
                        agentParams.params.lambda = rand([funcOptions.nStarts,1]);
                    case {'twoLinear', 'twoExponential'}
                        agentParams.params.lambda = rand([funcOptions.nStarts,1]);
                        agentParams.params.gamma = rand([funcOptions.nStarts,1]);
                end

            case 'epsilon-greedy'
                agentParams.params.epsilon= rand([funcOptions.nStarts,1]);

            case 'epsilon-softmax'
                agentParams.params.beta = exprnd(3,[funcOptions.nStarts,1]);
                agentParams.params.epsilon= rand([funcOptions.nStarts,1]);
        end

        % set up learning rates for average RR
        if modelOptions.learnRho || strcmp(modelOptions.rhoFunction, 'rw')
            agentParams.params.alphaRho = rand([funcOptions.nStarts,1]);
        end

        % set up learning rates for patch RR
        if modelOptions.learnPatchRR
            agentParams.params.alphaPatch = rand([funcOptions.nStarts,1]);
        end

        agentParams.params = struct2table(agentParams.params);
        agentParams.names = agentParams.params.Properties.VariableNames; 

        % set lower and upper bounds for fmincon search
        lb = array2table([0,0,0,0,-20,-20]); lb.Properties.VariableNames =  {'beta', 'epsilon', 'alphaRho', 'alphaPatch', 'lambda', 'gamma'};
        ub = array2table([50,50,1,1,20,20]); ub.Properties.VariableNames =  {'beta', 'epsilon', 'alphaRho', 'alphaPatch', 'lambda', 'gamma'};

        agentParams.lb = table2array(lb(:,agentParams.names)); % only get the parameters we need for this model
        agentParams.ub = table2array(ub(:,agentParams.names));
        agentParams.nParams = numel(agentParams.ub);
end
