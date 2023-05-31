function [out] = simulateM1_MVT(params, Block, Env)
%% Set up
PatchOrder = Env.PatchOrder{Block}; % specify which patches they see (high, medium, low quality) depending on environment

% reinforcement larning parameters
AlphaPatchRR = params(1); % Q learning rate
%AlphaRho = 0; % average RR learning rate

% initialise variables
Action = zeros(Env.BlockTime+1,1); % what action taken 
Reward = zeros(Env.BlockTime+1,1); % the reward earned in each state in the block
PatchRR = zeros(Env.BlockTime+1,1); % the reward rate in each state in the block 
PatchRPE = zeros(Env.BlockTime+1,1); % reward prediction error
RhoRPE = zeros(Env.BlockTime+1,1); % reward prediction error
%Rho = zeros(Env.BlockTime+1, 1); % estimated averageRR
Rho = ones(Env.BlockTime+1, 1)*Env.OptimalAverageRR{Block}; % initialise to actual averageRR
EstimatedPatchRR = zeros(Env.BlockTime+1, 1); % estimated patchRR

% initialise possible actions
Leave = 1;
Stay = 2;

NumObservations = 0; % how many updates to parameters (excluding leave states and arrival states) 

% set the first actions for the first patch
PatchNumber = 0; % start outside the patch
Arrive = 1; % start by arriving at new patch
Action(1) = Stay; % set first action to stay

% run model
for ii = 1:Env.BlockTime % for each second in the environment

    if Action(ii) == Stay % take action to stay

        if Arrive % if arriving in the patch
            T = 1; % time in patch
            PatchNumber = PatchNumber + 1; % patch number increases
            PatchType = PatchOrder(PatchNumber); % select patch type
            Arrive = 0; % no longer arriving
        end

        % observe reward in current state
        Reward(ii) = Env.R(T,PatchType); % reward depends on time in patch and patch type
        PatchRR(ii) = Reward(ii)/Env.TimeStep; % Reward Rate - this is the same as the reward, according to TimeStep. 
        
        PatchRPE(ii) = Reward(ii) - EstimatedPatchRR(ii);
        %RhoRPE(ii) = Reward(ii) - Rho(ii);

        % what is the estimated patch reward rate 
        EstimatedPatchRR(ii+1) = EstimatedPatchRR(ii) + AlphaPatchRR * PatchRPE(ii);
        %Rho(ii+1) = Rho(ii) + AlphaRho * RhoRPE(ii);

        Action(ii+1) = double(EstimatedPatchRR(ii+1) > Rho(ii+1)) + 1;    

        % if the next action is to leave 
        if Action(ii+1) == Leave
            t = 1; % if next action is leave, then reset travel time counter
            LeavingTime(PatchNumber,1) = T; % log the patch leaving time 
            LeavingRR(PatchNumber,1) = PatchRR(ii); % log the patch leaving reward rate 
        end

        T = T+Env.TimeStep; % time in patch increases
        NumObservations = NumObservations+1;

    elseif Action(ii) == Leave % take action to leave
        T = 0; % not in a patch anymore
        if t == Env.TravelTime % if on the last second of travelling
            Action(ii+1) = Stay;
            Arrive = 1; % about to arrive in new patch
        elseif t < Env.TravelTime % if still travelling
            Action(ii+1) = Leave;
        end

        % update 
        Reward(ii) = 0; % not getting anything during travel
        PatchRR(ii) = 0; % patch reward rate;
       
        PatchRPE(ii) = Reward(ii) - EstimatedPatchRR(ii);
        %RhoRPE(ii) = Reward(ii) - Rho(ii);

        % what is the estimated patch reward rate 
        EstimatedPatchRR(ii+1) = EstimatedPatchRR(ii) + AlphaPatchRR * PatchRPE(ii);
        %Rho(ii+1) = Rho(ii) + AlphaRho * RhoRPE(ii);

        t = t+1; % increase time spent travelling
    end
end

% store variables 
out.EstimatedPatchRR = EstimatedPatchRR(1:Env.BlockTime);
out.Rho = Rho(1:Env.BlockTime);
out.Action = Action(1:Env.BlockTime);
out.Reward = Reward(1:Env.BlockTime);
out.PatchRR = PatchRR(1:Env.BlockTime);
out.PatchRPE = PatchRPE(1:Env.BlockTime);
out.RhoRPE = RhoRPE(1:Env.BlockTime);
out.NumObservations = NumObservations;
out.LeavingTime = LeavingTime;
out.LeavingRR = LeavingRR;
out.PatchOrder = PatchOrder(1:length(LeavingTime));
end
