function [NegLogLikelihood, out] = NLL_MVT_softmax_dynamicbeta_timestep(params, Env, SubjData)

%% Transform their data into timesteps 
Env.TimeStep = params(2);

%% transform their data into timesteps 
% specify actions
Leave = 1;
Stay = 2;

LT = SubjData.summary.leaveT; % pull out leaving times - note that this will do it in the correct block order for the participant [1 2] or [2 1]
LT = round(LT/Env.TimeStep)*Env.TimeStep; % round to nearest timestep (state precision)

% transform leaving times into stay/leave actions for each state
for ii = 1:length(LT)
    a{ii} = repelem([Stay Leave], round([LT(ii) Env.TravelTime]./Env.TimeStep));
end
Action = cat(2, a{:})'; % concatenate all patches into one long block

PatchOrder = SubjData.summary.patch; % pull out patch order
BlockType = SubjData.summary.env; % Block type
BlockTime = size(Action, 1);

%% re-estimate reward based on timestep parameter
Env.RR = []; % refresh
Env.R = []; % refresh

for j=1:length(Env.InitialYield) % for each yield type
    fun = @(T) Env.InitialYield(j)*exp(-T*0.075); %Exponential function
    Env.RR(:,j) = fun(0:Env.TimeStep:max(LT)); % for max leave time for each participant
    k = 0:Env.TimeStep:max(LT); %this is set for a maximum of block_time secs in the patch, with loop time according to time_step
    for i = 1:length(k)
        Env.R(i,j) = integral(fun,k(i),k(i)+Env.TimeStep); % This creates matrix with each cell reward earned in time_step epoch at that reward rate
    end
end

%% Set up

Lambda = params(1); % weight for softmax temperature 

% initialise variables
PAction = zeros(BlockTime,2); % probability of selecting leave or stay
Reward = zeros(BlockTime,1); % the reward earned in each state in the block
PatchRR = zeros(BlockTime,1); % the reward rate in each state in the block 
experiencedAvgRR = zeros(BlockTime+1,1); % real experienced average RR

% initialise values - these will depend on model type
PAction(1, :) = [0 1]; % set first probabilities (starting in patch)

LogLikelihood = 0; 

% set the first actions for the first patch
PatchNumber = 0; % start outside the patch
Arrive = 1; % start by arriving at new patch

sumReward = 0; 
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
        PatchRR(ii) = Env.RR(N,PatchType); % Reward Rate - this is the same as the reward, according to TimeStep. 
        
        sumReward = Reward(ii)+sumReward;
        %experiencedAvgRR(ii) = sumReward/(ii);
        %experiencedAvgRR(ii) = Env.OptimalAverageRR{BlockType(PatchNumber)};
        experiencedAvgRR(ii) = SubjData.experiencedAvgRR(BlockType(PatchNumber));

        Beta = Lambda * 1/experiencedAvgRR(ii); % beta function based on avgRR

        PAction(ii+1,:) = CorrectedSoftmax([0, PatchRR(ii)], Beta);
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

        if Action(ii+1) == Stay % if entering new patch on next timstep
            PAction(ii+1,:) = [0 1]; % force staying on next action
            Arrive = 1; % about to arrive in new patch
        else 
            PAction(ii+1,:) = [1 0]; % force leaving for duration of travel time
        end

        % update 
        Reward(ii) = 0; % not getting anything during travel
        PatchRR(ii) = 0; % patch reward rate
        sumReward = Reward(ii)+sumReward;
        %experiencedAvgRR(ii) = sumReward/ii;
        %experiencedAvgRR(ii) = Env.OptimalAverageRR{BlockType(PatchNumber)};
        experiencedAvgRR(ii) = SubjData.experiencedAvgRR(BlockType(PatchNumber));

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
out.experiencedAvgRR = experiencedAvgRR; 
out.numObservations = sum(Action == 2); 

end
