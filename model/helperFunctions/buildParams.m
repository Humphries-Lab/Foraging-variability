function [agentParams] = buildParams(task,modelOptions,funcOptions)

switch funcOptions.type
    case 'simulate_new'

        if strcmp(funcOptions.blockPresentation,'combined')
            funcOptions.params.poor = funcOptions.params.rich; % only read in one set of parameters
        end

        % convert rich parameters to table to index by name 
        T = array2table(repmat(funcOptions.params.rich,[funcOptions.nSim,1])); 
        T.Properties.VariableNames = {'beta', 'epsilon', 'alphaRho', 'alphaPatch', 'lambda', 'gamma', 'bias'}; tmp{1} = T;
        
        % convert poor parameters to table to index by name 
        T = array2table(repmat(funcOptions.params.poor,[funcOptions.nSim,1])); 
        T.Properties.VariableNames = {'beta', 'epsilon', 'alphaRho', 'alphaPatch', 'lambda', 'gamma', 'bias'}; tmp{2} = T; clear T

        for iE = 1:task.nEnviron

            % set up action policy parameters
            switch modelOptions.actionPolicy

                case 'softmax'

                    switch modelOptions.betaFunction
                        case 'fit'
                            agentParams.params{iE}.beta = tmp{iE}.beta;
                        case {'scalar', 'exponential', 'hyperbolic'}
                            agentParams.params{iE}.lambda = tmp{iE}.lambda;
                        case {'twoScalar', 'twoExponential', 'twoHyperbolic'}
                            agentParams.params{iE}.lambda = tmp{iE}.lambda;
                            agentParams.params{iE}.gamma = tmp{iE}.gamma;
                    end

                case 'epsilon-greedy'
                    agentParams.params{iE}.epsilon = tmp{iE}.epsilon;

                case 'epsilon-softmax'
                    agentParams.params{iE}.beta = tmp{iE}.beta;
                    agentParams.params{iE}.epsilon = tmp{iE}.epsilon;
            end

            % set up learning rates for average RR
            if modelOptions.learnRho
                agentParams.params{iE}.alphaRho = tmp{iE}.alphaRho;
            end

            % set up learning rates for patch RR
            if modelOptions.learnPatchRR
                agentParams.params{iE}.alphaPatch = tmp{iE}.alphaPatch;
            end

            % set up bias (intercept) parameter for softmax
            if modelOptions.bias
                agentParams.params{iE}.bias = tmp{iE}.bias;
            end

            agentParams.params{iE} = struct2table(agentParams.params{iE});

        end
            agentParams.names = agentParams.params{1}.Properties.VariableNames; 

    case 'simulate_fit'

        if contains(funcOptions.blockPresentation, 'separate')
            load(sprintf('../data/fitting_data_240104/fitting_results_%s_M%d_%s', funcOptions.blockPresentation, modelOptions.modelNumber, funcOptions.study), 'minNLLFitParams_rich', 'minNLLFitParams_poor');
            agentParams.params{1} = minNLLFitParams_rich;
            agentParams.params{2} = minNLLFitParams_poor;
            agentParams.names = agentParams.params{1}.Properties.VariableNames;

        elseif contains(funcOptions.blockPresentation, 'combined')
            load(sprintf('../data/fitting_data_240104/fitting_results_%s_M%d_%s', funcOptions.blockPresentation, modelOptions.modelNumber, funcOptions.study), 'minNLLFitParams');
            agentParams.params{1} = minNLLFitParams;
            agentParams.params{2} = minNLLFitParams; % replicate for '2nd' simulation block
            agentParams.names = agentParams.params{1}.Properties.VariableNames;
        end

    case 'fit'
        % generate set of parameters to start fmincon search

        % set up action policy parameters
        switch modelOptions.actionPolicy

            case 'softmax'

                switch modelOptions.betaFunction
                    case 'fit'
                        agentParams.params.beta = exprnd(1,[funcOptions.nSim,1]);
                    case {'scalar', 'exponential', 'hyperbolic'}
                        agentParams.params.lambda = rand([funcOptions.nSim,1]);
                    case {'twoScalar', 'twoExponential', 'twoHyperbolic'}
                        agentParams.params.lambda = rand([funcOptions.nSim,1]);
                        agentParams.params.gamma = rand([funcOptions.nSim,1]);
                end

            case 'epsilon-greedy'
                agentParams.params.epsilon = rand([funcOptions.nSim,1]);

            case 'epsilon-softmax'
                agentParams.params.beta = exprnd(3,[funcOptions.nSim,1]);
                agentParams.params.epsilon= rand([funcOptions.nSim,1]);
        end

        % set up learning rates for average RR
        if modelOptions.learnRho
            agentParams.params.alphaRho = rand([funcOptions.nSim,1]);
        end

        % set up learning rates for patch RR
        if modelOptions.learnPatchRR
            agentParams.params.alphaPatch = rand([funcOptions.nSim,1]);
        end

                % set up bias (intercept) parameter
        if modelOptions.bias
            agentParams.params.bias = rand([funcOptions.nSim,1]);
        end

        agentParams.params = struct2table(agentParams.params);
        agentParams.names = agentParams.params.Properties.VariableNames; 

        % set lower and upper bounds for fmincon search
        lb = array2table([0,0,0,0,0,0, -20]); lb.Properties.VariableNames =  {'beta', 'epsilon', 'alphaRho', 'alphaPatch', 'lambda', 'gamma', 'bias'};
        ub = array2table([50,50,1,1,50,50, 20]); ub.Properties.VariableNames =  {'beta', 'epsilon', 'alphaRho', 'alphaPatch', 'lambda', 'gamma', 'bias'};

        agentParams.lb = table2array(lb(:,agentParams.names)); % only get the parameters we need for this model
        agentParams.ub = table2array(ub(:,agentParams.names));
        agentParams.nParams = numel(agentParams.ub);
end
