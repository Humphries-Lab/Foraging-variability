function [negLL,results,out] = simulate_MVT_model(task,model,agent,agentParams,funcOptions)

% Foraging variability project 
% Emma Scholey, 16 January 2025 

% simulate_MVT_model is main simulation script for all models

% INPUTS:
%   task: provides the task variables specific to each dataset
%
%   model: specifies which model we're simulating/fitting
%
%   agent: data container for each subject, in the format specified by the
%   buildData.m function. Includes real subject actions if fitting, or
%   empty if not simulating
%
%   agentParams: parameters specific to the model we're using
%
%   funcOptions: are we simulating or fitting? Required for telling the
%   script when to switch blocks. 

% OUTPUTS:
%   negLL: negative log likelihood for model fitting
%   results: summarised results: leaving time for each patch visited, for further plotting
%   out: fine-grained results: probabilities and actions taken at each timestep 

%% Set up

% convert parameter array to table to allow indexing
params = array2table(agentParams); params.Properties.VariableNames = model.paramNames;

% possible actions
leave = 1;
stay = 2;

% set up the subject's data containers
action = agent.data.action;
rho = agent.data.rho; % average reward rate
estPatchRR = agent.data.estPatchRR;

action(1) = stay; % start by staying in a patch
patchN = 0; % start outside the patch
arrive = 1; % start by arriving at new patch

logLikelihood = 0; % for fitting

for ii = 1:agent.nStates % for each second in the task

    if action(ii) == stay % take action to stay

        if arrive % if arriving in the patch

            timeNow = task.timeStep; % time in patch
            n = 1; % index (separate to time step)
            patchN = patchN + 1; % patch number increases
            patchType = agent.patchOrder(patchN); % select patch type
            arrive = 0; % no longer arriving

            if strcmp(funcOptions.type, 'fit') % if fitting data, need to index current block within the script (e.g. if we're fitting the whole session)
                currentEnv = agent.env(patchN);
            else
                currentEnv = agent.currentEnv;
            end

            if patchN == 1 || (strcmp(funcOptions.type, 'fit') && agent.switchIndex(patchN)==1) % if at first patch of a new block, check new environment and index new average reward rate

                switch model.startRho
                    case 'mvt'
                        rho(ii) = task.optAvgRR(currentEnv); % initialise rho to the MVT optimal avg reward rate for that environment
                    case 'block'
                        rho(ii) = agent.experiencedAvgRR(currentEnv); % set rho to agents' estimate of experiencedAvgRR.
                    case 'null'
                        rho(ii) = 0; % rho is not estimated, as fitting softmax temperature directly
                end
            end
        end

        switch task.rewardFunction
            case 'exponential'
                reward = reward_at_t_exp(timeNow,task.r0(patchType),task.decayRate); % reward depends on time in patch and patch type
            case 'linear'
                reward = reward_at_t_linear(timeNow,task.r0(patchType),task.decayRate); % reward depends on time in patch and patch type
        end

        % update MVT decision variables
        if model.learnPatchRR == 1
            patchRPE = reward - estPatchRR(ii); % estimate patch reward rate with delta learning rule
            estPatchRR(ii+1) = estPatchRR(ii) + params.alphaPatch * patchRPE;
        else
            estPatchRR(ii+1) = reward; % just assume full knowledge of the patch decay
        end

        if model.learnRho
            rhoRPE = reward - rho(ii); % estimate average reward rate with delta learning rule
            rho(ii+1) = rho(ii) + params.alphaRho * rhoRPE;
        else
            rho(ii+1) = rho(ii); % just assume full knowledge of the average reward rate
        end

        % compute beta
        switch model.betaFunction
            case 'single'
                current_beta = params.beta;
            case 'separate'
                betas = [params.beta_rich, params.beta_poor];
                current_beta = betas(currentEnv);
            case 'scalar' % remaining equations link average reward rate to beta 
                current_beta = params.lambda * 1/rho(ii+1);
            case 'exponential'
                current_beta = rho(ii+1)^params.gamma;
            case 'hyperbolic'
                current_beta = 1/(1 + params.lambda * rho(ii+1));
            case 'scalar_exp'
                current_beta = params.lambda * rho(ii+1)^params.gamma;
        end

        trackBeta(ii,1) = current_beta; % to plot trajectories after

        % set bias depending on environment and model
        switch model.bias
            case 'null'
                bias = 0; % set bias to 0 for models that don't require it
            case 'single'
                bias = params.bias;
            case 'separate'
                biases = [params.bias_rich, params.bias_poor];
                bias = biases(currentEnv);
        end

        % calculate action probabilities
        switch model.actionPolicy
            case 'softmax'
                pLeave = p_leave_softmax(estPatchRR(ii+1), current_beta, bias, 0); % don't pass rho to the function - use instead to modulate beta
        end
        
        trackP(ii,1) = pLeave;
        pAction = [pLeave, 1-pLeave];

        if strcmp(funcOptions.type, 'fit') % if fitting data
            pSelected = pAction(action(ii+1)); % what is probability of action they actually took
            pSelected(pSelected == 0) = eps(0); % prevent log(pselected = 0) going to infinity
            logLikelihood = logLikelihood + log(pSelected); % update log likelihood
        else % if simulating, then simulate their actions based on probabilities
            action(ii+1) = discreteinvrnd(pAction,1,1) ; 
        end

        % if the next action is to leave
        if action(ii+1) == leave
            t = task.timeStep; % if next action is leave, then reset travel time counter
            leaveT(patchN,1) = timeNow; % log the patch leaving time
            cE(patchN,1) = currentEnv; % log the environment for this patch (for recovery purposes)
        end

        timeNow = timeNow+task.timeStep; % time in patch increases
        n = n+1; % index increases

    elseif action(ii) == leave % take action to leave

        timeNow = 0; % not in a patch anymore

        if t == task.travelTime(currentEnv) % if on the last second of travelling
            action(ii+1) = stay;
            arrive = 1; % about to arrive in new patch
        elseif t < task.travelTime(currentEnv) % if still travelling
            action(ii+1) = leave;
        end

        reward = 0; % no reward during travelling

        % update MVT decision variables
        if model.learnPatchRR == 1
            patchRPE = reward - estPatchRR(ii);
            estPatchRR(ii+1) = estPatchRR(ii) + params.alphaPatch * patchRPE;
        else
            estPatchRR(ii+1) = reward; % just assume full knowledge of the patch decay
        end

        if model.learnRho
            rhoRPE = reward - rho(ii);
            rho(ii+1) = rho(ii) + params.alphaRho * rhoRPE;
        else
            rho(ii+1) = rho(ii);
        end

        % even though beta isn't used for anything here, track it anyway
        % compute beta
        switch model.betaFunction
            case 'single'
                current_beta = params.beta;
            case 'separate'
                betas = [params.beta_rich, params.beta_poor];
                current_beta = betas(currentEnv);
            case 'scalar'
                current_beta = params.lambda * 1/rho(ii+1);
            case 'exponential'
                current_beta = rho(ii+1)^params.gamma;
            case 'hyperbolic'
                current_beta = 1/(1 + params.lambda * rho(ii+1));
            case 'scalar_exp'
                current_beta = params.lambda * rho(ii+1)^params.gamma;
        end

        trackBeta(ii,1) = current_beta;
        t = t+task.timeStep; % increase time spent travelling
    end
end


out.rho = rho;
out.estPatchRR = estPatchRR;
out.action = action;
out.beta = trackBeta;
out.pLeave = trackP;

results.patchOrder = agent.patchOrder(1:numel(leaveT));
results.leaveT = leaveT;
results.env = cE; % add on extra row to account for recovery
results = struct2table(results);

negLL = -logLikelihood;

