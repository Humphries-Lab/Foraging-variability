function [agentParams] = buildParams(task,modelOptions,funcOptions)

switch funcOptions.type
    case 'simulate_new'

        % convert rich parameters to table to index by name
        T = array2table(repmat(funcOptions.params,[funcOptions.nSim,1]));
        T.Properties.VariableNames = {'beta','beta_rich', 'beta_poor','epsilon', 'alphaRho', 'alphaPatch', 'lambda', 'gamma', 'bias', 'bias_rich','bias_poor'};
        tmp = T;

        % set up action policy parameters
        switch modelOptions.actionPolicy

            case 'softmax'

                switch modelOptions.betaFunction
                    case 'single'
                        agentParams.params.beta = tmp.beta;
                    case 'separate'
                        agentParams.params.beta_rich = tmp.beta_rich;
                        agentParams.params.beta_poor = tmp.beta_poor; % only take single value of beta for both environments
                    case {'scalar', 'exponential', 'hyperbolic'}
                        agentParams.params.lambda = tmp.lambda;
                    case {'twoScalar', 'twoExponential', 'twoHyperbolic'}
                        agentParams.params.lambda = tmp.lambda;
                        agentParams.params.gamma = tmp.gamma;
                end

            case 'epsilon-greedy'
                agentParams.params.epsilon = tmp.epsilon;

            case 'epsilon-softmax'
                agentParams.params.beta = tmp.beta;
                agentParams.params.epsilon = tmp.epsilon;
        end

        % set up learning rates for average RR
        if modelOptions.learnRho
            agentParams.params.alphaRho = tmp.alphaRho;
        end

        % set up learning rates for patch RR
        if modelOptions.learnPatchRR
            agentParams.params.alphaPatch = tmp.alphaPatch;
        end

        % set up bias (intercept) parameter for softmax
        switch modelOptions.bias
            case 'single'
                agentParams.params.bias = tmp.bias;
            case 'separate'
                agentParams.params.bias_rich = tmp.bias_rich;
                agentParams.params.bias_poor = tmp.bias_poor;
        end

        agentParams.params = struct2table(agentParams.params);

        agentParams.names = agentParams.params.Properties.VariableNames;

    case 'simulate_fit'
        load(sprintf('../data/fitting_data_240117/fitting_results_M%d_%s', modelOptions.modelNumber, funcOptions.study), 'minNLLFitParams');
        agentParams.params = minNLLFitParams;
        agentParams.names = agentParams.params.Properties.VariableNames;

    case 'fit'
        % generate set of parameters to start fmincon search

        % set up action policy parameters
        switch modelOptions.actionPolicy

            case 'softmax'

                switch modelOptions.betaFunction
                    case 'single'
                        agentParams.params.beta = exprnd(1,[funcOptions.nSim,1]);
                    case 'separate'
                        agentParams.params.beta_rich = exprnd(1,[funcOptions.nSim,1]);
                        agentParams.params.beta_poor = exprnd(1,[funcOptions.nSim,1]);
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
        switch modelOptions.bias
            case 'single'
                agentParams.params.bias = rand([funcOptions.nSim,1]);
            case 'separate'
                agentParams.params.bias_rich = rand([funcOptions.nSim,1]);
                agentParams.params.bias_poor = rand([funcOptions.nSim,1]);
        end

        agentParams.params = struct2table(agentParams.params);
        agentParams.names = agentParams.params.Properties.VariableNames;

        % set lower and upper bounds for fmincon search
        lb = array2table([0,0,0,0,0,0,0,0,-20,-20,-20]); lb.Properties.VariableNames =  {'beta','beta_rich', 'beta_poor','epsilon', 'alphaRho', 'alphaPatch', 'lambda', 'gamma', 'bias', 'bias_rich','bias_poor'};
        ub = array2table([50,50,50,50,1,1,50,50,20,20,20]); ub.Properties.VariableNames =  {'beta', 'beta_rich', 'beta_poor','epsilon', 'alphaRho', 'alphaPatch', 'lambda', 'gamma', 'bias','bias_rich','bias_poor'};

        agentParams.lb = table2array(lb(:,agentParams.names)); % only get the parameters we need for this model
        agentParams.ub = table2array(ub(:,agentParams.names));
        agentParams.nParams = numel(agentParams.ub);
end
