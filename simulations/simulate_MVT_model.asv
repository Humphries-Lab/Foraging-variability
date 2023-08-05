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
newBlock = 0; % only required when fitting combined blocks together 

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
                if diff(agent.env(patchN-1:patchN)) % if fitting data, need to tell model when a new block is encountered 
                    newBlock = 1; 
                end
            else
                currentBlock = agent.currentBlock; 
            end
        end

        % observe reward in current state
        reward = reward_at_t_exp(timeNow,task.r0(patchType),task.decayRate); % reward depends on time in patch and patch type
        
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
                if timeNow == 1 && (patchN == 1 || newBlock == 1) % if in first patch, then initialise experiencedAvgRR in right range for the task
                    experiencedAvgRR = task.optAvgRR(currentBlock);
                end
                experiencedAvgRR = experiencedAvgRR + params.alphaRho*(reward - experiencedAvgRR); % rescorla wagner rule to update avgRR estimate
        end

        % compute beta 
        switch model.betaFunction
            case 'fit'
                beta = params.beta;
            case 'linear'
                beta = params.lambda * 1/experiencedAvgRR;
            case 'exponential'
                beta = 1/(experiencedAvgRR^params.lambda);
            case 'twoLinear'
                beta = params.lambda * 1/experiencedAvgRR*params.gamma;
            case 'twoExponential'
                beta = 1/(experiencedAvgRR^params.lambda * params.gamma);
        end


        % make choice
        switch model.actionPolicy
            case 'softmax'
                pAction = softmaxConstrain(estPatchRR(ii+1), rho(ii+1), beta);
            case 'epsilon-greedy'
                pAction(leave) = p_leave_egreedy(estPatchRR(ii+1), rho(ii+1), beta);
                pAction(stay) = 1 - pAction(leave);
            case 'epsilon-softmax'
                pAction = epsilon_softmaxConstrain(estPatchRR(ii+1), rho(ii+1), beta);
            case 'deterministic'
                pAction(stay) = estPatchRR(ii+1) > rho(ii+1); 
                pAction(leave) = 1 - pAction(stay);
        end

        if strcmp(funcOptions.type, 'fit') % if fitting data
            pSelected = pAction(action(ii+1)); % what is probability of action they actually took
            logLikelihood = logLikelihood + log(pSelected); % update log likelihood
        else
            action(ii+1) = discreteinvrnd(pAction,1,1) ; % simulate their actions based on probabilities
        end

        % if the next action is to leave
        if action(ii+1) == leave
            t = task.timeStep; % if next action is leave, then reset travel time counter
            leaveT(patchN,1) = timeNow; % log the patch leaving time
        end

        timeNow = timeNow+task.timeStep; % time in patch increases
        n = n+1; % index increases

    elseif action(ii) == leave % take action to leave

        timeNow = 0; % not in a patch anymore

        if t == task.travelTime % if on the last second of travelling
            action(ii+1) = stay;
            arrive = 1; % about to arrive in new patch
        elseif t < task.travelTime % if still travelling
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

results.patchOrder = agent.patchOrder(1:numel(leaveT));
results.leaveT = leaveT;
results = struct2table(results);

negLL = -logLikelihood;

