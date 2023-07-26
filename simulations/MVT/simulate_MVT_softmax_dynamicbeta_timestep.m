function [out] = simulate_MVT_softmax_dynamicbeta_timestep(params, Env, startValues)
%% Set up

Lambda = params(1); % softmax temperature
Env.TimeStep = params(2);

% initialise possible actions
Leave = 1;
Stay = 2;

% initialise variables
PAction = zeros(Env.BlockTime/Env.TimeStep+1,2); % probability of selecting leave or stay
Action = zeros(Env.BlockTime/Env.TimeStep+1,1); % what action taken
Reward = zeros(Env.BlockTime/Env.TimeStep+1,1); % the reward earned in each state in the block
PatchRR = zeros(Env.BlockTime/Env.TimeStep+1,1); % the reward rate in each state in the block
%experiencedAvgRR = zeros(Env.BlockTime/Env.TimeStep+1,1); % real experienced average RR
Beta = zeros(Env.BlockTime/Env.TimeStep+1,1);

%experiencedAvgRR = Env.OptimalAverageRR{Env.CurrentBlock}*ones(Env.BlockTime/Env.TimeStep+1,1);
experiencedAvgRR = Env.experiencedAvgRR(Env.CurrentBlock)*ones(Env.BlockTime/Env.TimeStep+1,1);

PatchOrder = Env.BlockPatchOrder;

% set the first actions for the first patch
PAction(1, :) = [0 1]; % set first probabilities (starting in patch)
PatchNumber = 0; % start outside the patch
Arrive = 1; % start by arriving at new patch
Action(1) = Stay; % set first action to stay

sumReward = startValues.sumReward;
timeInPreviousBlock = startValues.timeInPreviousBlock;

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

        sumReward = Reward(ii)+sumReward;
        %experiencedAvgRR(ii) = sumReward/(ii+timeInPreviousBlock);
        Beta(ii) = Lambda * 1/experiencedAvgRR(ii); % beta function based on avgRR

        PAction(ii+1,:) = CorrectedSoftmax([0, PatchRR(ii)], Beta(ii));

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
        if round(t,1) == Env.TravelTime % if on the last second of travelling
            PAction(ii+1,:) = [0 1]; % force staying on next action
            Action(ii+1) = Stay;
            Arrive = 1; % about to arrive in new patch
        elseif round(t,1) < Env.TravelTime % if still travelling
            PAction(ii+1,:) = [1 0]; % force leaving for duration of travel time
            Action(ii+1) = Leave;
        end

        % update
        Reward(ii) = 0; % not getting anything during travel
        PatchRR(ii) = 0; % patch reward rate
        sumReward = Reward(ii)+sumReward;
        %experiencedAvgRR(ii) = sumReward/(ii+timeInPreviousBlock);
        Beta(ii) = Lambda * 1/experiencedAvgRR(ii); % beta function based on avgRR

        t = t+Env.TimeStep; % increase time spent travelling
    end
end

% store variables
out.PAction = PAction(1:Env.BlockTime/Env.TimeStep,:);
out.Action = Action(1:Env.BlockTime/Env.TimeStep);
out.Reward = Reward(1:Env.BlockTime/Env.TimeStep);
out.PatchRR = PatchRR(1:Env.BlockTime/Env.TimeStep);
out.LeavingTime = LeavingTime;
out.LeavingRR = LeavingRR;
out.PatchOrder = PatchOrder(1:length(LeavingTime));
out.sumReward = sumReward;
out.experiencedAvgRR = experiencedAvgRR(1:Env.BlockTime/Env.TimeStep);
out.Beta = Beta(1:Env.BlockTime/Env.TimeStep);

out.Rho = zeros(Env.BlockTime+1, 1);
out.EstimatedPatchRR = zeros(Env.BlockTime+1, 1); % we include this as part of the output data to account for startValues

end
