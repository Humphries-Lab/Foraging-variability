
function [out] = Simulate_RLOn_Model(RLparams,BlockOrder,Task,Model)

% This function simulates ON-policy RL model.
% 
%   INPUTS:
%   RLparams: RL parameters
%           -Alpha -R learning rate
%           -Beta -softmax temperature
%           -AlphaRho -average RR learning rate
%           -Epsilon (in case of e-greedy model)
%           
%   BlockOrder: Block order presentation ([1 2] rich-poor or [2 1]
%   poor-rich)
%
%   Task: a structure contaning Task parameters
%           -TimeStep in seconds
%           -TravelTime: travel time between patches in seconds
%           -BlockNum: number of blocks
%           -BlockTime: duration of each block in seconds
%           -Duration: task duration (depending on combined session or block presented separately)
%           -PatchNum: number of patch types
%           -PatchReward: The reward offered in each state of a patch type
%           -Block_POrder: Patch order per block (environment type)
%           -MTVOptimal: MVT optimal leaving time prediction for each block in seconds 
%
%   Model: a structure containing information about the type of RL model
%           -Num: Model number as specified in excel document
%           -Patch_sequence: indicates whether Rstay sequence is multiple(per patch type) or single (regardless of patch type)
%           -Leave_sequence: indicates whether Rleave has states (one per travel time) or is stateless (single value updated in each travel time)
%           -Next_patch: prediccion of next patch (only if multiple patch
%           sequence). Type of predictions considered are 'previous'(same as patch just left), 
%           'Godlike'(next patch is known),'Stochastic'(guess randomly or based on the most frequent patch encountered in the block)
%           -Action_selection: indicates whether next action is selected based on e-greedy or softmax (stochastic)function.
%           -Block_presentation: indicates whether to run both blocks separately or a combined session.
%   
%   OUTPUT:
%   out:    Structure with the results of the RL model simulation. Most important parameters are:
%           -LeavingTime: the leaving time resulting from the simulation
%           -LeavingRR: the leaving reward rate at the time of leaving each
%           patch
%           -Rstay and Rleave values
%           -NextPatch guess if multiple patch sequence
%           -NumObservations
%           -R,Rho,RPE,Reward,PatchRR, PAction, Action.
%
%   Version: 17/07/2023
%   Francisca Perea Perez
%% SET UP VARIABLES
% Possible actions
Leave = 1;
Stay = 2;

% RL parameters
AlphaR = RLparams(1);   % R learning rate
AlphaRho = RLparams(2); % average RR learning rate
switch Model.Action_selection   %if Stochastic selection I need Beta, if e-greedy I need Epsilon.
    case 'Stochastic'
        Beta = RLparams(3);   % softmax temperature
    otherwise
        Epsilon= RLparams(3); % epsilon
        n=10;                       % if e-greedy choice, we force the model to stay at least n seconds in the first patch(es) when Rstay=0
end
    
% Initialise variables (one cell per block) 
R = repmat({zeros(Task.BlockTime+1,2)},1,2);         % [R(leave), R(stay)]  
PAction = repmat({zeros(Task.BlockTime+1,2)},1,2);   % probability of selecting leave or stay (per block)
Action = repmat({zeros(Task.BlockTime+1,1)},1,2);    % what action taken
Reward = repmat({zeros(Task.BlockTime+1,1)},1,2);    % the reward earned in each state in the block
PatchRR = repmat({zeros(Task.BlockTime+1,1)},1,2);   % the reward rate in each state in the block
RPE = repmat({zeros(Task.BlockTime+1,1)},1,2);       % reward prediction error
Rho = repmat({zeros(Task.BlockTime+1,1)},1,2);       % estimated averageRR

% Rstay and Rleave - Based on Model type
switch Model.Patch_sequence
    case 'Multiple'
        RStay = zeros(Task.Duration+1,Task.PatchNum);       % R per patch type
        NextPatch = cell([1 2]);                            % Store my predictions of next patch while leaving
        PPatch= ones(1,Task.PatchNum)*(1/Task.PatchNum);    % Prior probability of each patch type
        NPatchMemory= 10;                                   % Number of patches to consider if short memory
    case 'Single'
        RStay = zeros(Task.Duration+1, 1);                    % single R over all patches
end
switch Model.Leave_sequence
    case 'Separate'
        RLeave = zeros(Task.TravelTime+1,1);                 % all seconds in leave states are represented
    case 'Stateless'
        RLeave = 0;                                                % single leave state
end

% Reward normalisation
maxR = max(max(Task.PatchReward)); % Maximum possible reward

% Info about Parameters updates
NumObservations = repmat({0},1,2); % how many updates excluding leave states and arrival states

% Set up storage for results (one cell per block)
LeavingTime = cell([1 2]);      % Per block
LeavingRR = cell ([1 2]);       % Per block

% Block Presentation
PatchesSeen= [];                                % Patches seen throughout the block or the combined session
if strcmpi(Model.Block_presentation,'Separate') %if separate simulation, store Rstay and RLeave per block
    RStayBlock= cell([1 2]);
    RLeaveBlock= cell([1 2]);
end


%% Run either main combined session or separate blocks

for b=1:Task.BlockNum         % for each block

    % Update new block environment & first actions
    BlockType = BlockOrder(b);                    % Next block
    PatchOrder = Task.Block_POrder(:,BlockType);   % set new patch order
    PatchNumber = 0;                                     % start outside the patch
    Arrive = 1;                                          % start by arriving at new patch
    Action{BlockType}(1)= Stay;                          % set first action to stay
    PAction{BlockType}(1,:) = [0 1];                     % 100% probability of staying
    
    for s = 1:Task.BlockTime     % for each state (second) of the block   
                      
        %IF STAYING
        if Action{BlockType}(s) == Stay 

            if Arrive % if arriving in the patch
                T = 1;                                          % time in patch
                PatchNumber = PatchNumber + 1;                  % patch number increases
                PatchType = PatchOrder(PatchNumber);            % select patch type
                Arrive = 0;                                     % no longer arriving
                PatchesSeen = cat(2,PatchesSeen,PatchType);     % update patches seen 
                if strcmpi(Model.Patch_sequence,'Multiple') 
                   R{BlockType}(s,Stay) = RStay(T, PatchType);  % change value of first state according to the patch I am in        
                end                
            end

        % Observe reward in current state
        Reward{BlockType}(s) = Task.PatchReward(T,PatchType);    % reward depends on time in patch and patch type
        PatchRR{BlockType}(s) = Reward{BlockType}(s)/Task.TimeStep;         % Reward Rate - Reward per timestep.

        % Set R values for next state to select action
            switch Model.Patch_sequence
                case 'Multiple'
                    R{BlockType}(s+1,:) = [RLeave(1), RStay(T+1,PatchType)]; % select Rstay based on patch_type
                case 'Single'
                    R{BlockType}(s+1,:) = [RLeave(1), RStay(T+1)]; % no patch type to select
            end
            
            
        % Next action selection
           switch Model.Action_selection   
                case 'Stochastic'
                    PAction{BlockType}(s+1,:) = exp(R{BlockType}(s+1,:) * Beta) ./ sum(exp(Beta*R{BlockType}(s+1,:))); % probability of next action
                    Action{BlockType}(s+1) = discreteinvrnd(PAction{BlockType}(s+1,:),1,1);                 % Choose next action depending on probability PAction                    
                otherwise     % Epsilon-greedy choice
                   if rand< Epsilon
                       Action{BlockType}(s+1) = discreteinvrnd([0.5 0.5], 1, 1);  % random value
                   else
                       if R{BlockType}(s+1,Stay)==0 && T<=n   % set the greedy action to stay for n seconds at the beggining (first patch(es)) 
                          % when all Rstay values are 0 to avoid the model from getting stuck in leaving action                       
                           Action{BlockType}(s+1) = Stay;
                       else                           
                           if R{BlockType}(s+1,1)~= R{BlockType}(s+1,2)
                                Action{BlockType}(s+1) = find(R{BlockType}(s+1,:)== max(R{BlockType}(s+1,:)));    % Choose the maximum value action
                           else % if both values are the same, flip a coin
                                Action{BlockType}(s+1) = discreteinvrnd([0.5 0.5], 1, 1);  % random value
                           end                           
                       end
                   end                   
                   
           end
                              
        % if Next action is leave
            if Action{BlockType}(s+1) == Leave
                t = 1; % if next action is leave, then reset travel time counter
                LeavingTime{BlockType}(PatchNumber) = T; % log the patch leaving time
                LeavingRR{BlockType}(PatchNumber) = PatchRR{BlockType}(s); % log the patch leaving reward rate
                
                % what does agent think the NEXT patch is?                
                if exist('NextPatch','var')
                    switch Model.Next_patch
                        case 'Previous'
                            NextPatch{BlockType}(PatchNumber) = PatchType;                 % assumes next patch is the same as the current one
                        case 'Godlike'
                            NextPatch{BlockType}(PatchNumber) = PatchOrder(PatchNumber+1); % knows what the next patch is
                        
                        case 'Stochastic'                            
                            if strcmpi(Model.Memory,'Short')                                
                                % Calculate the probability of each patch type based on last patches "NPatchMemory" experienced
                                PPatch= sum(PatchesSeen(max(1,size(PatchesSeen,2)-(NPatchMemory-1)):end)'== (1:Task.PatchNum))./ min(size(PatchesSeen,2),NPatchMemory);    
                                NextPatch{BlockType}(PatchNumber) = discreteinvrnd(PPatch,1,1);               % Choose Next patch based on their probability 
                                
                            elseif strcmpi(Model.Memory,'Block')
                                PPatch= sum(PatchOrder(1:PatchNumber)== (1:Task.PatchNum))./  PatchNumber;    % Calculate the probability of each patch type based on experience in current block
                                NextPatch{BlockType}(PatchNumber) = discreteinvrnd(PPatch,1,1);               % Choose Next patch based on their probability 
                                
                            else % if memory covers the entire session (either combined or separate blocks)
                                PPatch= sum(PatchesSeen(1:end)'== (1:Task.PatchNum))./ size(PatchesSeen,2);   % Calculate the probability of each patch type based on experience in block or entire session
                                NextPatch{BlockType}(PatchNumber) = discreteinvrnd(PPatch,1,1);               % Choose Next patch based on their probability                                 
                            end              
                            
                        otherwise % ?-greedy approach
                            if (rand < 0.1)  % for 10% of the time, choose randomly
                                NextPatch{BlockType}(PatchNumber) = randi([1 3]);
                            else % else assume most frequent patch                                
                                if strcmpi(Model.Memory,'Short')
                                    NextPatch{BlockType}(PatchNumber) = mode(PatchesSeen(max(1,size(PatchesSeen,2)-(NPatchMemory-1)):end));                                
                                elseif strcmpi(Model.Memory,'Block')
                                    NextPatch{BlockType}(PatchNumber) = mode(PatchOrder(1:PatchNumber)); % in current block                                
                                else % if memory covers the entire session (either combined or separate blocks)
                                    NextPatch{BlockType}(PatchNumber) = mode(PatchesSeen(1:end));
                                end
                            end
                                                                
                         end
                        
                   end
                end
            
            
        % Update Rstay          
        RPE{BlockType}(s) = Reward{BlockType}(s)/maxR - Rho{BlockType}(s) + R{BlockType}(s+1, Action{BlockType}(s+1)) - R{BlockType}(s, Action{BlockType}(s));
        Rho{BlockType}(s+1) = Rho{BlockType}(s) + AlphaRho * RPE{BlockType}(s); % update estimate of average RR
        switch Model.Patch_sequence
            case 'Multiple'
                RStay(T,PatchType) = RStay(T,PatchType) + AlphaR * RPE{BlockType}(s);
            case 'Single'
                RStay(T) = RStay(T) + AlphaR * RPE{BlockType}(s);
        end
       
        % Time in patch and N.of observations increase
        T = T+Task.TimeStep; 
        NumObservations{BlockType} = NumObservations{BlockType}+1;
        
        
        % IF LEAVING
        elseif Action{BlockType}(s) == Leave 
            T = 0; % not in a patch anymore            
            
            % Set R values for next state
            switch Model.Patch_sequence
                case 'Multiple'
                    if strcmpi(Model.Leave_sequence,'Separate')
                        R{BlockType}(s+1,:) = [RLeave(t+1), RStay(1,NextPatch{BlockType}(PatchNumber))]; %Rleave separate states & Qstay multiple sequence
                    else
                        R{BlockType}(s+1,:) = [RLeave, RStay(1,NextPatch{BlockType}(PatchNumber))];    % Rleave Stateless & Rstay multiple sequence                       
                    end
                case 'Single'
                    if strcmpi(Model.Leave_sequence,'Separate')
                            R{BlockType}(s+1,:) = [RLeave(t+1), RStay(1)];    % Rleave separate states & Rstay single sequence
                    else
                            R{BlockType}(s+1,:) = [RLeave, RStay(1)];         % Rleave Stateless & Rstay single sequence
                    end
            end
            
            % Next action selection depending on travel time                   
            if t == Task.TravelTime % if on the last second of travelling
                PAction{BlockType}(s+1,:) = [0 1]; % force staying on next action
                Action{BlockType}(s+1) = Stay;
                Arrive = 1; % about to arrive in new patch
            elseif t < Task.TravelTime % if still travelling
                PAction{BlockType}(s+1,:) = [1 0]; % force leaving for duration of travel time
                Action{BlockType}(s+1) = Leave;    

            end
            
            % Update RLeave
            Reward{BlockType}(s) = 0; % not getting any reward during travel
            PatchRR{BlockType}(s) = 0; % patch reward rate;
            RPE{BlockType}(s) = (Reward{BlockType}(s)/maxR - Rho{BlockType}(s)) + R{BlockType}(s+1, Action{BlockType}(s+1)) - R{BlockType}(s, Action{BlockType}(s)); % calculate RPE for this state
            Rho{BlockType}(s+1) = Rho{BlockType}(s) + AlphaRho * RPE{BlockType}(s); % update estimate of average RR
            switch Model.Leave_sequence
                case 'Separate'
                    RLeave(t) = RLeave(t) + AlphaR * RPE{BlockType}(s); % update Rleave based on individual states
                case'Stateless'
                    RLeave = RLeave + AlphaR * RPE{BlockType}(s); % update single Rleave value
            end
            
            % Increase time spent travelling
            t = t+1; 
            
        end
%           if s==300
%               fprintf('chequeo')
%           end
   
     end

        if  (b < Task.BlockNum) && (strcmpi(Model.Block_presentation,'Combined')) % If Combined simulation
            
            % last block's state values are passed to the first state of the next block
            R{BlockOrder(b+1)}(1,:)= R{BlockOrder(b)}(end,:);
            Rho{BlockOrder(b+1)}(1)= Rho{BlockOrder(b)}(end);
            %NumObservations keeps incrementing in the new block
            NumObservations{BlockOrder(b+1)}= NumObservations{BlockOrder(b)}; 
        
        elseif (b < Task.BlockNum) && (strcmpi(Model.Block_presentation,'Separate')) % If Separate simulation
            
            %Save Rstay and Rleave values from first block
            RStayBlock{BlockOrder(b)}= RStay(1:max(sum(RStay~=0)),:);
            RLeaveBlock{BlockOrder(b)}= RLeave;            
            % Reset Patches seen before running the next block
            PatchesSeen= [];            
            % Reset Rstay and Rleave values before running next block                    
            switch Model.Patch_sequence
                case 'Multiple'
                    RStay = zeros(Task.Duration+1,Task.PatchNum);   % RStay per patch type
                case 'Single'
                    RStay = zeros(Task.Duration+1, 1);              % single R over all patches
            end
            switch Model.Leave_sequence
                case 'Separate'
                    RLeave = zeros(Task.TravelTime+1,1);            % all seconds in leave states are represented
                case 'Stateless'
                    RLeave = 0;                                     % single leave state
            end                    
                    
        end                 
            
end

%% Store variables
% Save RStay and RLeave
switch Model.Block_presentation
    case 'Combined'
        out.RStay = RStay(1:max(sum(RStay~=0)),:); % limit R_stay tables - remove values not used   
        out.RLeave = RLeave;
        out.NumObservations = NumObservations{BlockOrder(b)};
    case 'Separate'  
        %Store the last block 
        RStayBlock{BlockOrder(b)}= RStay(1:max(sum(RStay~=0)),:); 
        RLeaveBlock{BlockOrder(b)}= RLeave;
        %Save
        out.RStay= RStayBlock;
        out.RLeave= RLeaveBlock;
        out.NumObservations = NumObservations;
end

if strcmpi(Model.Patch_sequence,'Multiple')
    out.NextPatch = NextPatch;
end

out.R = R; 
out.Rho = cellfun(@(x) x*maxR,Rho,'un',false);
out.PAction = PAction;
out.Action = Action;
out.Reward = Reward;
out.PatchRR = PatchRR;
out.RPE = RPE;
out.LeavingTime = LeavingTime;
out.LeavingRR = LeavingRR;
out.BlockOrder = BlockOrder; 
out.PatchOrder = [{Task.Block_POrder((1:size(LeavingTime{1},2)+1),1)}, {Task.Block_POrder((1:size(LeavingTime{2},2)+1),2)}];% log extra patch for each environment (+1), just to allow for the last patch they are in
out.PresentationType= Model.Block_presentation;     
      
                                 
    end
            
    

    
    


