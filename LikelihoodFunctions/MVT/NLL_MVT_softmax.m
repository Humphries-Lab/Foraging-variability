function [NegLogLikelihood, out] = NLL_MVT_softmax(params, Env, SubjData)
%% Set up
PatchOrder = SubjData.PatchOrder; % specify which patches they see (high, medium, low quality) depending on environment
Action = SubjData.Action;

BlockTime = size(Action, 1);

Beta = params(1); % softmax temperature

% initialise variables
PAction = zeros(BlockTime,2); % probability of selecting leave or stay
Reward = zeros(BlockTime,1); % the reward earned in each state in the block
PatchRR = zeros(BlockTime,1); % the reward rate in each state in the block 

% initialise values - these will depend on model type
PAction(1, :) = [0 1]; % set first probabilities (starting in patch)

% initialise possible actions
Leave = 1;
Stay = 2;

LogLikelihood = 0; 

% set the first actions for the first patch
PatchNumber = 0; % start outside the patch
Arrive = 1; % start by arriving at new patch

% run model
for ii = 1:BlockTime-1 % for each subject action

    if Action(ii) == Stay % take action to stay

        if Arrive % if arriving in the patch
            T = Env.TimeStep; % time in patch
            N = 1;
            PatchNumber = PatchNumber + 1; % patch number increases
            PatchType = PatchOrder(PatchNumber); % select patch type
            Arrive = 0; % no longer arriving
        end

        % observe reward in current state
        Reward(ii) = Env.R(N,PatchType); % reward depends on time in patch and patch type
        PatchRR(ii) = Reward(ii)/Env.TimeStep; % Reward Rate - this is the same as the reward, according to TimeStep. 
        
        %PAction(ii+1,:) = CorrectedSoftmax([0, PatchRR(ii)], Beta);
        PAction(ii+1,Leave) = p_leave_softmax(PatchRR(ii), Beta); 
        PAction(ii+1,Stay) = 1 - PAction(ii+1,Leave);
        
        pSelected = PAction(ii+1, Action(ii+1)); % what did they actually do next, and what PAction does the model estimate
        LogLikelihood = LogLikelihood + log(pSelected); % update log likelihood

        % if the next action is to leave 
        if Action(ii+1) == Leave
            t = Env.TimeStep; % if next action is leave, then reset travel time counter
            LeavingTime(PatchNumber,1) = T; % log the patch leaving time
            LeavingRR(PatchNumber,1) = PatchRR(ii); % log the patch leaving reward rate
        end

        T = T+Env.TimeStep; % time in patch increases
        N = N+1;

    elseif Action(ii) == Leave % take action to leave
        T = 0; 

        if round(t,1) == Env.TravelTime % if on the last timestep of travelling. Issue with floating point integer requires round
            PAction(ii+1,:) = [0 1]; % force staying on next action
            Arrive = 1; % about to arrive in new patch
        elseif round(t,1) < Env.TravelTime % if still travelling
            PAction(ii+1,:) = [1 0]; % force leaving for duration of travel time
        end

        % update 
        Reward(ii) = 0; % not getting anything during travel
        PatchRR(ii) = 0; % patch reward rate

        t = t+Env.TimeStep; % increase time spent travelling
    end
end

NegLogLikelihood = -LogLikelihood;

% store variables 
out.PAction = PAction(1:BlockTime,:);
out.Action = Action(1:BlockTime);
out.Reward = Reward(1:BlockTime);
out.PatchRR = PatchRR(1:BlockTime-1);
out.LeavingTime = LeavingTime;
out.LeavingRR = LeavingRR;
out.PatchOrder = PatchOrder(1:length(LeavingTime));
end
