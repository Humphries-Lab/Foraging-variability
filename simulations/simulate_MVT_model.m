function [negLogLikelihood,out,summary] = simulate_MVT_model(task,model,agent,agentParams)

%% Set up

% possible actions
leave = 1;
stay = 2;

out = agent.data; 

% set the first actions for the first patch
out.pAction(1, :) = [0 1]; % set first probabilities (starting in patch)
out.action(1) = stay; % set first action to stay

patchN = 0; % start outside the patch
arrive = 1; % start by arriving at new patch

patchOrder = agent.blockPatchOrder;

negLoglikelihood = 0;

for ii = 1:numel(out.action) % for each second in the task

    if out.action(ii) == stay % take action to stay

        if arrive % if arriving in the patch
            timeNow = task.timeStep; % time in patch
            n = 1; % index (separate to time step)
            patchN = patchN + 1; % patch number increases
            patchType = patchOrder(patchN); % select patch type
            arrive = 0; % no longer arriving
        end

        % observe reward in current state
        r = reward_at_t_exp(timeNow,task.r0(patchType),task.decayRate); % reward depends on time in patch and patch type
        
        out.reward(ii) = r; 
        out.patchRR(ii) = out.reward(ii)/task.timeStep;

        patchRPE = r - out.estPatchRR(ii);
        rhoRPE = r - out.rho(ii);

        % what is the estimated patch reward rate
        out.estPatchRR(ii+1) = out.estPatchRR(ii) + agentParams.alphaPatch * patchRPE;
        out.rho(ii+1) = out.rho(ii) + agentParams.alphaRho * rhoRPE;

        out.pAction(ii+1,:) = CorrectedSoftmax([out.rho(ii+1), out.estPatchRR(ii+1)], agentParams.beta);

        % choose next action: discreteinvrnd will output 1 (leave) or 2 (stay), depending on probability distribution of PAction (1-p, p).
        out.action(ii+1) = discreteinvrnd(out.pAction(ii+1,:),1,1) ;

        % if the next action is to leave
        if out.action(ii+1) == leave
            t = task.timeStep; % if next action is leave, then reset travel time counter
            leaveT(patchN,1) = timeNow; % log the patch leaving time
        end

        timeNow = timeNow+task.timeStep; % time in patch increases
        n = n+1; % index increases
        %negLogLikelihood = negLogLikelihood + log(pSelected);

    elseif out.action(ii) == leave % take action to leave
        timeNow = 0; % not in a patch anymore
        if t == task.travelTime % if on the last second of travelling
            out.pAction(ii+1,:) = [0 1]; % force staying on next action
            out.action(ii+1) = stay;
            arrive = 1; % about to arrive in new patch
        elseif t < task.travelTime % if still travelling
            out.pAction(ii+1,:) = [1 0]; % force leaving for duration of travel time
            out.action(ii+1) = leave;
        end

        % update
        r = 0; 
        out.reward(ii) = r; % not getting anything during travel
        out.patchRR(ii) = r/task.timeStep; % patch reward rate;

        patchRPE = r - out.estPatchRR(ii);
        rhoRPE = r - out.rho(ii);

        % what is the estimated patch reward rate
        out.estPatchRR(ii+1) = out.estPatchRR(ii) + agentParams.alphaPatch * patchRPE;
        out.rho(ii+1) = out.rho(ii) + agentParams.alphaRho * rhoRPE;

        t = t+task.timeStep; % increase time spent travelling
    end

end

summary.patchOrder = patchOrder(1:numel(leaveT))';
summary.leaveT = leaveT; 
summary.block = repelem(agent.blockType,numel(leaveT))';
summary = struct2table(summary);
end
