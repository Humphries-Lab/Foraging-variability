function [negLL,results,out] = simulate_MVT_model(task,model,agent,agentParams,funcOptions)

%% Set up

% convert parameter array to table to allow indexing
params = array2table(agentParams); params.Properties.VariableNames = model.paramNames;

% possible actions
leave = 1;
stay = 2;

action = agent.data.action;
rho = agent.data.rho;
estPatchRR = agent.data.estPatchRR;

action(1) = stay; % start by staying in a patch
patchN = 0; % start outside the patch
arrive = 1; % start by arriving at new patch

logLikelihood = 0;

for ii = 1:agent.nStates % for each second in the task

    if action(ii) == stay % take action to stay

        if arrive % if arriving in the patch

            timeNow = task.timeStep; % time in patch
            n = 1; % index (separate to time step)
            patchN = patchN + 1; % patch number increases
            patchType = agent.patchOrder(patchN); % select patch type
            arrive = 0; % no longer arriving

            if strcmp(funcOptions.type, 'fit') % if fitting data, need to index current block within the script (e.g. if we're fitting the whole session)
                currentBlock = agent.env(patchN);
            else
                currentBlock = agent.currentBlock;
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
            patchRPE = reward - estPatchRR(ii);
            estPatchRR(ii+1) = estPatchRR(ii) + params.alphaPatch * patchRPE;
        else
            estPatchRR(ii+1) = reward; % just assume full knowledge of the patch decay
        end

        if model.learnRho == 1
            rhoRPE = reward - rho(ii);
            rho(ii+1) = rho(ii) + params.alphaRho * rhoRPE;
        else
            rho(ii+1) = rho(ii);
        end

        % compute experiencedAvgRR for beta function models
        switch model.rhoFunction
            case 'mvt'
                experiencedAvgRR = task.optAvgRR(currentBlock); % MVT optimal
            case 'block'
                experiencedAvgRR = agent.experiencedAvgRR(currentBlock); % real experienced avgRR
            case 'rw'
                if timeNow == 1 % if first started the session, then initialise experiencedAvgRR in right range for the task
                    experiencedAvgRR = task.optAvgRR(currentBlock);
                end
                experiencedAvgRR = experiencedAvgRR + params.alphaRho*(reward - experiencedAvgRR); % rescorla wagner rule to update avgRR estimate
        end

        % make choice
        switch model.actionPolicy
            case 'softmax'
                % compute beta
                switch model.betaFunction
                    case 'fit'
                        beta = params.beta;
                    case 'scalar'
                        beta = params.lambda * 1/experiencedAvgRR;
                    case 'exponential'
                        beta = 1/(experiencedAvgRR^params.lambda);
                    case 'hyperbolic'
                        beta = 1/(1 + params.lambda * experiencedAvgRR);
                    case 'twoScalar'
                        beta = params.lambda * 1/(params.gamma * experiencedAvgRR);
                    case 'twoExponential'
                        beta = params.lambda * 1/(experiencedAvgRR^params.gamma);
                    case 'twoHyperbolic'
                        beta = params.lambda * 1/(1 + params.lambda * experiencedAvgRR);
                end

                pAction = softmaxConstrain(estPatchRR(ii+1), rho(ii+1), beta);

            case 'epsilon-greedy'
                pAction = [params.epsilon, 1-params.epsilon]; % is this correct when rho does not equal 0? i.e. not just a coin flip on the patch reward rate now
            case 'epsilon-softmax'
                pAction = epsilon_softmaxConstrain(estPatchRR(ii+1), rho(ii+1), params.beta, params.epsilon);
            case 'deterministic'
                pAction(stay) = double(estPatchRR(ii+1) > rho(ii+1));
                pAction(leave) = 1 - pAction(stay);
        end

        if strcmp(funcOptions.type, 'fit') % if fitting data
            pSelected = pAction(action(ii+1)); % what is probability of action they actually took
            pSelected(pSelected == 0) = eps(0); % prevent log(pselected = 0) going to infinity 
            logLikelihood = logLikelihood + log(pSelected); % update log likelihood
        else
            action(ii+1) = discreteinvrnd(pAction,1,1) ; % simulate their actions based on probabilities
        end

        % if the next action is to leave
        if action(ii+1) == leave
            t = task.timeStep; % if next action is leave, then reset travel time counter
            leaveT(patchN,1) = timeNow; % log the patch leaving time
            cB(patchN,1) = currentBlock; % log the environment for this patch (for recovery purposes) 
        end

        timeNow = timeNow+task.timeStep; % time in patch increases
        n = n+1; % index increases

    elseif action(ii) == leave % take action to leave

        timeNow = 0; % not in a patch anymore

        if t == task.travelTime(currentBlock) % if on the last second of travelling
            action(ii+1) = stay;
            arrive = 1; % about to arrive in new patch
        elseif t < task.travelTime(currentBlock) % if still travelling
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

        if model.learnRho == 1
            rhoRPE = reward - rho(ii);
            rho(ii+1) = rho(ii) + params.alphaRho * rhoRPE;
        else
            rho(ii+1) = rho(ii);
        end

        t = t+task.timeStep; % increase time spent travelling
    end

end

out.rho = rho;
out.estPatchRR = estPatchRR;
out.action = action;

results.patchOrder = agent.patchOrder(1:numel(leaveT)+1); % account for extra patch they may have been in before task simulation ended
results.leaveT = leaveT;
results.env = [cB;currentBlock];

negLL = -logLikelihood;

