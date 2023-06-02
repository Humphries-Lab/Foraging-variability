function [out] = simulate_MVT_patchRR_RW_softmax(params, Env, startValues)
%% Set up

% reinforcement learning parameters
AlphaPatchRR = params(1);
Beta = params(2); % softmax temperature

% initialise possible actions
Leave = 1;
Stay = 2;

% initialise variables
PAction = zeros(Env.BlockTime/Env.TimeStep+1,2); % probability of selecting leave or stay
Action = zeros(Env.BlockTime/Env.TimeStep+1,1); % what action taken 
Reward = zeros(Env.BlockTime/Env.TimeStep+1,1); % the reward earned in each state in the block
PatchRR = zeros(Env.BlockTime/Env.TimeStep+1,1); % the reward rate in each state in the block 
EstimatedPatchRR = zeros(Env.BlockTime+1, 1); % estimated patchRR

EstimatedPatchRR(1) = startValues.EstimatedPatchRR;

PatchOrder = Env.BlockPatchOrder;

% set the first actions for the first patch
PAction(1, :) = [0 1]; % set first probabilities (starting in patch)
PatchNumber = 0; % start outside the patch
Arrive = 1; % start by arriving at new patch
Action(1) = Stay; % set first action to stay

% run model
for ii = 1:Env.BlockTime/Env.TimeStep % for each second in the environment

    if Action(ii) == Stay % take action to stay

        if Arrive % if arriving in the patch
            T = Env.TimeStep; % time in patch
            N = 1; % index (separate to time step)
            PatchNumber = PatchNumber + 1; % patch number increases
            PatchType = PatchOrder(PatchNumber); % select patch type
            Arrive = 0; % no longer arriving
        end

        % observe reward in current state
        Reward(ii) = Env.R(N,PatchType); % reward depends on time in patch and patch type
        PatchRR(ii) = Reward(ii)/Env.TimeStep;
        
        PatchRPE = Reward(ii) - EstimatedPatchRR(ii);

        % what is the estimated patch reward rate 
        EstimatedPatchRR(ii+1) = EstimatedPatchRR(ii) + AlphaPatchRR * PatchRPE;

        PAction(ii+1,:) = CorrectedSoftmax([0, EstimatedPatchRR(ii+1)], Beta);

        % choose next action: discreteinvrnd will output 1 (leave) or 2 (stay), depending on probability distribution of PAction (1-p, p). 
        Action(ii+1) = discreteinvrnd(PAction(ii+1,:),1,1) ;  
        
        % if the next action is to leave 
        if Action(ii+1) == Leave
            t = Env.TimeStep; % if next action is leave, then reset travel time counter
            LeavingTime(PatchNumber,1) = T; % log the patch leaving time 
            LeavingRR(PatchNumber,1) = PatchRR(ii); % log the patch leaving reward rate 
        end

        T = T+Env.TimeStep; % time in patch increases
        N = N+1; % index increases

    elseif Action(ii) == Leave % take action to leave
        T = 0; % not in a patch anymore
        if t == Env.TravelTime % if on the last second of travelling
            PAction(ii+1,:) = [0 1]; % force staying on next action
            Action(ii+1) = Stay;
            Arrive = 1; % about to arrive in new patch
        elseif t < Env.TravelTime % if still travelling
            PAction(ii+1,:) = [1 0]; % force leaving for duration of travel time
            Action(ii+1) = Leave;
        end

        % update 
        Reward(ii) = 0; % not getting anything during travel
        PatchRR(ii) = 0; % patch reward rate;
       
        PatchRPE = Reward(ii) - EstimatedPatchRR(ii);

        % what is the estimated patch reward rate 
        EstimatedPatchRR(ii+1) = EstimatedPatchRR(ii) + AlphaPatchRR * PatchRPE;

        t = t+Env.TimeStep; % increase time spent travelling
    end
end

% store variables 
out.EstimatedPatchRR = EstimatedPatchRR(1:Env.BlockTime);
out.PAction = PAction(1:Env.BlockTime/Env.TimeStep,:);
out.Action = Action(1:Env.BlockTime/Env.TimeStep);
out.Reward = Reward(1:Env.BlockTime/Env.TimeStep);
out.PatchRR = PatchRR(1:Env.BlockTime/Env.TimeStep);
out.LeavingTime = LeavingTime;
out.LeavingRR = LeavingRR;
out.PatchOrder = PatchOrder(1:length(LeavingTime));

out.Rho = zeros(Env.BlockTime+1, 1);
end
