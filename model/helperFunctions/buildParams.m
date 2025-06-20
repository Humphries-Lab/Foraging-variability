function [agentParams] = buildParams(modelOptions,funcOptions)

% Foraging variability project
% Emma Scholey, 16 January 2025

% BUILDPARAMS generates a set of model parameters for the stochastic
% foraging models.
% INPUTS:
%   modelOptions: based on the options specified in the foragingModelTable,
%   generate only the parameters required: reward-dependent beta (1 or 2), alpha (to learn
%   average reward rate or patch reward rate), lambda and gamma (to
%   translate average reward rate into gamma, and reward-independent bias
%   (1 or 2)
%
%   funcOptions: based on whether you are simulating with new
%   user-specified parameters (simulate_new), with pre-fit best model
%   parameters for each subject (simulate_fit), or fitting parameters to
%   the data and require a set of random initial starting seeds within
%   parameter bounds (fit)

switch funcOptions.type
    case 'simulate_new'

        T = array2table(repmat(funcOptions.params,[funcOptions.nSim,1]));
        T.Properties.VariableNames = {'beta','beta_rich', 'beta_poor','alphaRho', 'alphaPatch', 'lambda', 'gamma', 'bias', 'bias_rich','bias_poor'};
        tmp = T;

        switch modelOptions.betaFunction
            case 'single'
                agentParams.params.beta = tmp.beta;
            case 'separate'
                agentParams.params.beta_rich = tmp.beta_rich;
                agentParams.params.beta_poor = tmp.beta_poor; % only take single value of beta for both environments
            case {'scalar', 'hyperbolic'}
                agentParams.params.lambda = tmp.lambda;
            case 'exponential'
                agentParams.params.gamma = tmp.gamma;
            case 'scalar_exp'
                agentParams.params.lambda = tmp.lambda;
                agentParams.params.gamma = tmp.gamma;
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
        load(sprintf('../data/fitting_data/fitting_results_M%d_%s', modelOptions.modelNumber, funcOptions.study), 'minNLLFitParams');
        agentParams.params = minNLLFitParams;
        agentParams.names = agentParams.params.Properties.VariableNames;

    case 'fit'
        % generate set of parameters to start fmincon search

        % set up action policy parameters

        switch modelOptions.betaFunction
            case 'single'
                agentParams.params.beta = exprnd(1,[funcOptions.nSim,1]);
            case 'separate'
                agentParams.params.beta_rich = exprnd(1,[funcOptions.nSim,1]);
                agentParams.params.beta_poor = exprnd(1,[funcOptions.nSim,1]);
            case {'scalar', 'hyperbolic'}
                agentParams.params.lambda = 4*rand([funcOptions.nSim,1]);
            case 'exponential'
                agentParams.params.gamma = rand([funcOptions.nSim,1]);
            case 'scalar_exp'
                agentParams.params.lambda = 4*rand([funcOptions.nSim,1]);
                agentParams.params.gamma = rand([funcOptions.nSim,1]);
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
                agentParams.params.bias = -1*rand([funcOptions.nSim,1]);
            case 'separate'
                agentParams.params.bias_rich = -1*rand([funcOptions.nSim,1]);
                agentParams.params.bias_poor = -1*rand([funcOptions.nSim,1]);
        end

        agentParams.params = struct2table(agentParams.params);
        agentParams.names = agentParams.params.Properties.VariableNames;

        % set lower and upper bounds for fmincon search
        lb = array2table([0,0,0,0,0,0,-5,-60,-60,-60]); lb.Properties.VariableNames =  {'beta','beta_rich', 'beta_poor', 'alphaRho', 'alphaPatch', 'lambda', 'gamma', 'bias', 'bias_rich','bias_poor'};
        ub = array2table([50,50,50,1,1,inf,5,10,10,10]); ub.Properties.VariableNames =  {'beta', 'beta_rich', 'beta_poor', 'alphaRho', 'alphaPatch', 'lambda', 'gamma', 'bias','bias_rich','bias_poor'};

        agentParams.lb = table2array(lb(:,agentParams.names)); % only get the parameters we need for this model
        agentParams.ub = table2array(ub(:,agentParams.names));
        agentParams.nParams = numel(agentParams.ub);

    case 'recover'

        % load real participant data to fit pdf and generate simulated parameters
        load(sprintf('../data/fitting_data/fitting_results_M%d_%s', modelOptions.modelNumber, funcOptions.study), 'minNLLFitParams');

        % set up action policy parameters

        switch modelOptions.betaFunction
            case 'single'
               agentParams.params.beta = rand(funcOptions.nSim,1)*(mean(minNLLFitParams.beta) + 3*std(minNLLFitParams.beta));
                %agentParams.params.beta = rand(funcOptions.nSim,1)*max(minNLLFitParams.beta);

            case 'separate'
                agentParams.params.beta_rich = rand(funcOptions.nSim,1)*(mean(minNLLFitParams.beta_rich) + 3*std(minNLLFitParams.beta_rich));
                agentParams.params.beta_poor = rand(funcOptions.nSim,1)*(mean(minNLLFitParams.beta_poor) + 3*std(minNLLFitParams.beta_poor));
                %agentParams.params.beta_rich = rand(funcOptions.nSim,1)*max(minNLLFitParams.beta_rich);
                %agentParams.params.beta_poor = rand(funcOptions.nSim,1)*max(minNLLFitParams.beta_poor);
        end

        % set up bias (intercept) parameter
        switch modelOptions.bias
            case 'single'

                a = mean(minNLLFitParams.bias) - 3*std(minNLLFitParams.bias);
                b = mean(minNLLFitParams.bias) + 3*std(minNLLFitParams.bias);
                
                agentParams.params.bias = a+(b-a)*rand(funcOptions.nSim,1);
                %agentParams.params.bias = rand(funcOptions.nSim,1)* min(minNLLFitParams.bias);

            case 'separate'
                %agentParams.params.bias_rich = rand(funcOptions.nSim,1)*(mean(minNLLFitParams.bias_rich) - 3*std(minNLLFitParams.bias_rich));
                %agentParams.params.bias_poor = rand(funcOptions.nSim,1)*(mean(minNLLFitParams.bias_poor) - 3*std(minNLLFitParams.bias_poor));
                %agentParams.params.bias_rich = rand(funcOptions.nSim,1)* min(minNLLFitParams.bias_rich);
                %agentParams.params.bias_poor = rand(funcOptions.nSim,1)* min(minNLLFitParams.bias_poor);

                a = mean(minNLLFitParams.bias_rich) - 3*std(minNLLFitParams.bias_rich);
                b = mean(minNLLFitParams.bias_rich) + 3*std(minNLLFitParams.bias_rich);

                agentParams.params.bias_rich = a+(b-a)*rand(funcOptions.nSim,1);

                a = mean(minNLLFitParams.bias_poor) - 3*std(minNLLFitParams.bias_poor);
                b = mean(minNLLFitParams.bias_poor) + 3*std(minNLLFitParams.bias_poor);

                agentParams.params.bias_poor = a+(b-a)*rand(funcOptions.nSim,1);

        end

        agentParams.params = struct2table(agentParams.params);
        agentParams.names = agentParams.params.Properties.VariableNames;
        agentParams.nParams = numel(agentParams.names);


end
