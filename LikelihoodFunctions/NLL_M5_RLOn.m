function [NegLogLikelihood, out] = NLL_M5_RLOn(params, Env, SubjData)
%% Set up
PatchOrder = SubjData.PatchOrder; % specify which patches they see (high, medium, low quality) depending on environment
Action = SubjData.Action;

BlockTime = size(Action, 1);

% reinforcement learning parameters
AlphaQ = params(1); % Q learning rate
AlphaRho = params(2); % average RR learning rate
Beta = params(3); % softmax temperature

% initialise variables
Q = zeros(BlockTime,2); % [Q(leave), Q(stay)]
PAction = zeros(BlockTime,2); % probability of selecting leave or stay
Reward = zeros(BlockTime,1); % the reward earned in each state in the block
PatchRR = zeros(BlockTime,1); % the reward rate in each state in the block 
RPE = zeros(BlockTime,1); % reward prediction error
Rho = zeros(BlockTime, 1); % estimated averageRR
if SubjData.NextPatch == 0 % if no existing Next Patch predictions (i.e. not using simulated data) 
    NextPatch = zeros(BlockTime,1); % start from scratch 
else % otherwise, use the simulated data to ensure we correctly recover next patch predictions 
    NextPatch = SubjData.NextPatch; % take existing next patch predictions from simulated data
end

% initialise values - these will depend on model type
PAction(1, :) = [0 1]; % set first probabilities (starting in patch)
QStay = zeros(BlockTime,3); % Q-table for staying - Q per patch type - 'multiple' patch sequences
QLeave = zeros(Env.TravelTime+1,1); % Q-table for leaving - all seconds in leave states are represented - 'separate' leave sequence

% initialise possible actions
Leave = 1;
Stay = 2;


LogLikelihood = 0; 

% set the first actions for the first patch
PatchNumber = 0; % start outside the patch
Arrive = 1; % start by arriving at new patch

maxR = max(max(Env.R)); % what's the max possible reward, for normalisation 

% run model
for ii = 1:BlockTime-1 % for each second in the environment

    if Action(ii) == Stay % take action to stay

        if Arrive % if arriving in the patch
            T = 1; % time in patch
            PatchNumber = PatchNumber + 1; % patch number increases
            PatchType = PatchOrder(PatchNumber); % select patch type
            Arrive = 0; % no longer arriving
            Q(ii,Stay) = QStay(T, PatchType); % change value of first state based on what the patch actually is
        end

        % observe reward in current state
        Reward(ii) = Env.R(T,PatchType); % reward depends on time in patch and patch type
        PatchRR(ii) = Reward(ii)/Env.TimeStep; % Reward Rate - this is the same as the reward, according to TimeStep. 

        % what is the next state? 
        Q(ii+1,:) = [QLeave(1), QStay(T+1,PatchType)]; % s'

        PAction(ii+1,Stay) = softmaxStay(Beta, Q(ii+1, Stay), Q(ii+1,Leave)); % function to calculate PAction based on softmax
        PAction(ii+1,Leave) = 1 - PAction(ii+1,Stay); % p(Leave) is just inverse of p(Stay) 
        pSelected = PAction(ii+1, Action(ii+1)); % what did they actually do next, and what PAction does the model estimate
        % get rid of non-finite values for pSelected before updating likelihood
        if pSelected == 0
            pSelected = eps(0);
        end
        LogLikelihood = LogLikelihood + log(pSelected); % update log likelihood

        % if the next action is to leave 
        if Action(ii+1) == Leave
            t = 1; % if next action is leave, then reset travel time counter

            % what does agent think the NEXT patch is? 
            NextPatch(PatchNumber) = PatchType; % assume next patch will be same as previous patch

            LeavingTime(PatchNumber,1) = T; % log the patch leaving time
            LeavingRR(PatchNumber,1) = PatchRR(ii); % log the patch leaving reward rate
        end

        RPE(ii) = Reward(ii)/maxR - Rho(ii) + Q(ii+1, Action(ii+1)) - Q(ii, Action(ii));

        % update estimate of average RR
        Rho(ii+1) = Rho(ii) + AlphaRho * RPE(ii);

        % update Q table for staying
        QStay(T,PatchType) = QStay(T,PatchType) + AlphaQ * RPE(ii);

        T = T+Env.TimeStep; % time in patch increases

    elseif Action(ii) == Leave % take action to leave
        T = 0; 
        Q(ii+1,:) = [QLeave(t+1), QStay(1,NextPatch(PatchNumber))]; % what are Q value's for next state? 

        if t == Env.TravelTime % if on the last second of travelling
            PAction(ii+1,:) = [0 1]; % force staying on next action
            Arrive = 1; % about to arrive in new patch
        elseif t < Env.TravelTime % if still travelling
            PAction(ii+1,:) = [1 0]; % force leaving for duration of travel time
        end

        % update 
        Reward(ii) = 0; % not getting anything during travel
        PatchRR(ii) = 0; % patch reward rate;
        RPE(ii) = (Reward(ii)/maxR - Rho(ii)) + Q(ii+1, Action(ii+1)) - Q(ii, Action(ii)); % calculate RPE for this state
        Rho(ii+1) = Rho(ii) + AlphaRho * RPE(ii); % update estimate of average RR
        QLeave(t) = QLeave(t) + AlphaQ * RPE(ii); % update Q-leave table based on individual states

        t = t+1; % increase time spent travelling
    end
end

NegLogLikelihood = -LogLikelihood;

% store variables 
out.QStay = QStay(1:max(sum(QStay~=0)), :); % limit Q_stay tables - find longest column
out.QLeave = QLeave(1:Env.TravelTime);
out.Rho = Rho(1:BlockTime);
out.PAction = PAction(1:BlockTime,:);
out.Action = Action(1:BlockTime);
out.Q = Q(1:BlockTime,:);
out.Reward = Reward(1:BlockTime);
out.PatchRR = PatchRR(1:BlockTime);
out.RPE = RPE(1:BlockTime);
out.LeavingTime = LeavingTime;
out.LeavingRR = LeavingRR;
out.PatchOrder = PatchOrder(1:length(LeavingTime));
out.NextPatch = NextPatch(1:length(LeavingTime));
end
