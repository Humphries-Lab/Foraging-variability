function [out] = simulateM25_RLOn(params, Block, Env)
%% Set up
PatchOrder = Env.PatchOrder{Block}; % specify which patches they see (high, medium, low quality) depending on environment

% reinforcement learning parameters
AlphaQ = params(1); % Q learning rate
AlphaRho = params(2); % average RR learning rate
Beta = params(3); % softmax temperature

% initialise variables
Q = zeros(Env.BlockTime+1,2); % [Q(leave), Q(stay)]
PAction = zeros(Env.BlockTime+1,2); % probability of selecting leave or stay
Action = zeros(Env.BlockTime+1,1); % what action taken 
Reward = zeros(Env.BlockTime+1,1); % the reward earned in each state in the block
PatchRR = zeros(Env.BlockTime+1,1); % the reward rate in each state in the block 
RPE = zeros(Env.BlockTime+1,1); % reward prediction error
Rho = zeros(Env.BlockTime+1, 1); % estimated averageRR

% initialise values - these will depend on model type
PAction(1, :) = [0 1]; % set first probabilities (starting in patch)
QStay = zeros(Env.BlockTime+1,3); % Q-table for staying - Q per patch type - 'multiple' patch sequences
QLeave = zeros(Env.TravelTime+1,1); % Q-table for leaving - all seconds in leave states are represented - 'separate' leave sequence

% initialise possible actions
Leave = 1;
Stay = 2;

NumObservations = 0; % how many updates to parameters (excluding leave states and arrival states) 

% set the first actions for the first patch
PatchNumber = 0; % start outside the patch
Arrive = 1; % start by arriving at new patch
Action(1) = Stay; % set first action to stay

maxR = max(max(Env.R)); % what's the max possible reward, for normalisation 

% run model
for ii = 1:Env.BlockTime % for each second in the environment

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

        PAction(ii + 1,:) = exp(Q(ii+1, :) * Beta) ./ sum(exp(Beta*Q(ii+1, :))); % probability of next action 
         % choose next action: discreteinvrnd will output 1 (leave) or 2 (stay), depending on probability distribution of PAction (1-p, p). 
        Action(ii+1) = discreteinvrnd(PAction(ii+1,:),1,1) ;    

        % if the next action is to leave 
        if Action(ii+1) == Leave
            t = 1; % if next action is leave, then reset travel time counter
            % what does agent think the NEXT patch is? - stochastic patch prediction 
            if rand < 0.1 % for some % of the time, choose randomly
                NextPatch(PatchNumber) = randi([1 3]);
            else % else assume most frequent patch based on those previously experienced
                NextPatch(PatchNumber) = mode(PatchOrder(1:PatchNumber));
            end

            LeavingTime(PatchNumber,1) = T; % log the patch leaving time 
            LeavingRR(PatchNumber,1) = PatchRR(ii); % log the patch leaving reward rate 
        end

        RPE(ii) = Reward(ii)/maxR - Rho(ii) + Q(ii+1, Action(ii+1)) - Q(ii, Action(ii));

        % update estimate of average RR
        Rho(ii+1) = Rho(ii) + AlphaRho * RPE(ii);

        % update Q table for staying
        QStay(T,PatchType) = QStay(T,PatchType) + AlphaQ * RPE(ii);

        T = T+Env.TimeStep; % time in patch increases
        NumObservations = NumObservations+1;

    elseif Action(ii) == Leave % take action to leave
        T = 0; % not in a patch anymore
        Q(ii+1,:) = [QLeave(t+1), QStay(1,NextPatch(PatchNumber))]; % what are Q value's for next state? 

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
        RPE(ii) = (Reward(ii)/maxR - Rho(ii)) + Q(ii+1, Action(ii+1)) - Q(ii, Action(ii)); % calculate RPE for this state
        Rho(ii+1) = Rho(ii) + AlphaRho * RPE(ii); % update estimate of average RR
        QLeave(t) = QLeave(t) + AlphaQ * RPE(ii); % update Q-leave table based on individual states

        t = t+1; % increase time spent travelling
    end
end

% store variables 
out.QStay = QStay(1:max(sum(QStay~=0)), :); % limit Q_stay tables - find longest column
out.QLeave = QLeave(1:Env.TravelTime);
out.Rho = Rho(1:Env.BlockTime)*maxR;
out.PAction = PAction(1:Env.BlockTime,:);
out.Action = Action(1:Env.BlockTime);
out.Q = Q(1:Env.BlockTime,:);
out.Reward = Reward(1:Env.BlockTime);
out.PatchRR = PatchRR(1:Env.BlockTime);
out.RPE = RPE(1:Env.BlockTime);
out.NumObservations = NumObservations;
out.LeavingTime = LeavingTime;
out.LeavingRR = LeavingRR;
out.PatchOrder = PatchOrder(1:length(LeavingTime));
out.NextPatch = NextPatch;
end
