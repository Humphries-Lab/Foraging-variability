
function [out] = Simulate_RLOn_Model(RLparams,BlockOrder,Task,Model)

% This function simulates ON-policy RL model.
%
%   INPUTS:
%   RLparams: RL parameters
%           -Alpha -Q learning rate
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
%           -Patch_sequence: indicates whether Qstay sequence is multiple(per patch type) or single (regardless of patch type)
%           -Leave_sequence: indicates whether Qleave has states (one per travel time) or is stateless (single value updated in each travel time)
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
%           -Qstay and Qleave values
%           -NextPatch guess if multiple patch sequence
%           -NumObservations
%           -Q,Rho,RPE,Reward,PatchRR, PAction, Action.
%
%   Version: 17/07/2023
%   Francisca Perea Perez
%% SET UP VARIABLES
% Possible actions
Leave = 1;
Stay = 2;

% RL parameters
AlphaQ = RLparams(1);   % Q learning rate
AlphaRho = RLparams(2); % average RR learning rate
switch Model.Action_selection   %if Stochastic selection I need Beta, if e-greedy I need Epsilon.
    case 'Stochastic'
        Beta = RLparams(3);   % softmax temperature
    otherwise
        Epsilon= RLparams(3); % epsilon
        n=10;                       % if e-greedy choice, we force the model to stay at least n seconds in the first patch when Qstay=0
end

% Initialise variables (one cell per block)
Q = repmat({zeros(Task.BlockTime+1,2)},1,2);         % [Q(leave), Q(stay)]
PAction = repmat({zeros(Task.BlockTime+1,2)},1,2);   % probability of selecting leave or stay (per block)
Action = repmat({zeros(Task.BlockTime+1,1)},1,2);    % what action taken
Reward = repmat({zeros(Task.BlockTime+1,1)},1,2);    % the reward earned in each state in the block
PatchRR = repmat({zeros(Task.BlockTime+1,1)},1,2);   % the reward rate in each state in the block
RPE = repmat({zeros(Task.BlockTime+1,1)},1,2);       % reward prediction error
Rho = repmat({zeros(Task.BlockTime+1,1)},1,2);       % estimated averageRR

% Qstay and Qleave - Based on Model type
switch Model.Patch_sequence
    case 'Multiple'
        QStay = zeros(Task.Duration+1,Task.PatchNum);   % Q per patch type
        NextPatch = cell([1 2]);                                    % Store my predictions of next patch while leaving
    case 'Single'
        QStay = zeros(Task.Duration+1, 1);                    % single Q over all patches
end
switch Model.Leave_sequence
    case 'Separate'
        QLeave = zeros(Task.TravelTime+1,1);                 % all seconds in leave states are represented
    case 'Stateless'
        QLeave = 0;                                                % single leave state
end

% Reward normalisation
maxR = max(max(Task.PatchReward)); % Maximum possible reward

% Info about Parameters updates
NumObservations = repmat({0},1,2); % how many updates excluding leave states and arrival states

% Set up storage for results (one cell per block)
LeavingTime = cell([1 2]);      % Per block
LeavingRR = cell ([1 2]);       % Per block

% Block Presentation
if strcmpi(Model.Block_presentation,'Combined')
    PatchesSeen= [];                            % Patches seen throughout the combined session
else
    %if separate simulation, store Qstay and QLeave per block
    QStayBlock= cell([1 2]);
    QLeaveBlock= cell([1 2]);
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
                T = 1; % time in patch
                PatchNumber = PatchNumber + 1; % patch number increases
                PatchType = PatchOrder(PatchNumber); % select patch type
                Arrive = 0; % no longer arriving
                if strcmpi(Model.Patch_sequence,'Multiple')
                    Q{BlockType}(s,Stay) = QStay(T, PatchType);           % Change value of first state according to the patch I am in
                    if exist('PatchesSeen','var')
                        PatchesSeen = cat(2,PatchesSeen,PatchType); % Update patches seen in the entire session
                    end
                end
            end

            % Observe reward in current state
            Reward{BlockType}(s) = Task.PatchReward(T,PatchType);    % reward depends on time in patch and patch type
            PatchRR{BlockType}(s) = Reward{BlockType}(s)/Task.TimeStep;         % Reward Rate - Reward per timestep.

            % Set Q values for next state to select action
            switch Model.Patch_sequence
                case 'Multiple'
                    Q{BlockType}(s+1,:) = [QLeave(1), QStay(T+1,PatchType)]; % select Q_stay based on patch_type
                case 'Single'
                    Q{BlockType}(s+1,:) = [QLeave(1), QStay(T+1)]; % no patch type to select
            end


            % Next action selection
            switch Model.Action_selection
                case 'Stochastic'
                    PAction{BlockType}(s+1,:) = exp(Q{BlockType}(s+1,:) * Beta) ./ sum(exp(Beta*Q{BlockType}(s+1,:))); % probability of next action
                    Action{BlockType}(s+1) = discreteinvrnd(PAction{BlockType}(s+1,:),1,1);                 % Choose next action depending on probability PAction
                otherwise     % Epsilon-greedy choice
                    if rand< Epsilon
                        Action{BlockType}(s+1) = discreteinvrnd([0.5 0.5], 1, 1);  % random value
                    else
                        if Q{BlockType}(s+1,Stay)==0 && T<=n   % set the greedy action to stay for n seconds at the beggining (first patch(es))
                            % when all Qstay values are 0 to avoid the model from getting stuck in leaving action
                            Action{BlockType}(s+1) = Stay;
                        else
                            if Q{BlockType}(s+1,1)~= Q{BlockType}(s+1,2)
                                Action{BlockType}(s+1) = find(Q{BlockType}(s+1,:)== max(Q{BlockType}(s+1,:)));    % Choose the maximum value action
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
                            if rand < 0.1 % for some % of the time, choose randomly
                                NextPatch{BlockType}(PatchNumber) = randi([1 3]);
                            else % else assume most frequent patch
                                NextPatch{BlockType}(PatchNumber) = mode(PatchOrder(1:PatchNumber)); % in current block
                                %Other possible criteria
                                %                                 NextPatch{BlockType}(PatchNumber) = mode(PatchOrder(max(1,PatchNumber -10):PatchNumber)); % short memory: only considers last ten patches in block
                                %                                 if strcmpi(Model.Block_presentation,'Combined')
                                %                                     NextPatch{BlockType}(PatchNumber) = mode(PatchesSeen(1:end)) % most frequent patch in entire session
                                %                                     NextPatch{BlockType}(PatchNumber) = mode(PatchesSeen(max(1,size(PatchesSeen,2)-10):end)); % short memory: only considers last ten patches in session
                                %                                 end
                            end
                    end
                end
            end

            % Update Qstay
            RPE{BlockType}(s) = Reward{BlockType}(s)/maxR - Rho{BlockType}(s) + Q{BlockType}(s+1, Action{BlockType}(s+1)) - Q{BlockType}(s, Action{BlockType}(s));
            Rho{BlockType}(s+1) = Rho{BlockType}(s) + AlphaRho * RPE{BlockType}(s); % update estimate of average RR
            switch Model.Patch_sequence
                case 'Multiple'
                    QStay(T,PatchType) = QStay(T,PatchType) + AlphaQ * RPE{BlockType}(s);
                case 'Single'
                    QStay(T) = QStay(T) + AlphaQ * RPE{BlockType}(s);
            end

            % Time in patch and N.of observations increase
            T = T+Task.TimeStep;
            NumObservations{BlockType} = NumObservations{BlockType}+1;


            % IF LEAVING
        elseif Action{BlockType}(s) == Leave
            T = 0; % not in a patch anymore

            % Set Q values for next state
            switch Model.Patch_sequence
                case 'Multiple'
                    if strcmpi(Model.Leave_sequence,'Separate')
                        Q{BlockType}(s+1,:) = [QLeave(t+1), QStay(1,NextPatch{BlockType}(PatchNumber))]; %Qleave separate states & Qstay multiple sequence
                    else
                        Q{BlockType}(s+1,:) = [QLeave, QStay(1,NextPatch{BlockType}(PatchNumber))];    % Qleave Stateless & Qstay multiple sequence
                    end
                case 'Single'
                    if strcmpi(Model.Leave_sequence,'Separate')
                        Q{BlockType}(s+1,:) = [QLeave(t+1), QStay(1)];    % Qleave separate states & Qstay single sequence
                    else
                        Q{BlockType}(s+1,:) = [QLeave, QStay(1)];         % Qleave Stateless & Qstay single sequence
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

            % Update QLeave
            Reward{BlockType}(s) = 0; % not getting any reward during travel
            PatchRR{BlockType}(s) = 0; % patch reward rate;
            RPE{BlockType}(s) = (Reward{BlockType}(s)/maxR - Rho{BlockType}(s)) + Q{BlockType}(s+1, Action{BlockType}(s+1)) - Q{BlockType}(s, Action{BlockType}(s)); % calculate RPE for this state
            Rho{BlockType}(s+1) = Rho{BlockType}(s) + AlphaRho * RPE{BlockType}(s); % update estimate of average RR
            switch Model.Leave_sequence
                case 'Separate'
                    QLeave(t) = QLeave(t) + AlphaQ * RPE{BlockType}(s); % update Q-leave based on individual states
                case'Stateless'
                    QLeave = QLeave + AlphaQ * RPE{BlockType}(s); % update single Q_leave value
            end

            % Increase time spent travelling
            t = t+1;

        end

        %         if s==600 || s==1 || s==2
        %             sprintf('chequear');
        %         end

    end

    if  (b < Task.BlockNum) && (strcmpi(Model.Block_presentation,'Combined')) % If Combined simulation

        % last block's state values are passed to the first state of the next block
        Q{BlockOrder(b+1)}(1,:)= Q{BlockOrder(b)}(end,:);
        Rho{BlockOrder(b+1)}(1)= Rho{BlockOrder(b)}(end);
        %NumObservations keeps incrementing in the new block
        NumObservations{BlockOrder(b+1)}= NumObservations{BlockOrder(b)};

    elseif (b < Task.BlockNum) && (strcmpi(Model.Block_presentation,'Separate')) % If Separate simulation

        %Save Qstay and Qleave values from first block
        QStayBlock{BlockOrder(b)}= QStay(1:max(sum(QStay~=0)),:);
        QLeaveBlock{BlockOrder(b)}= QLeave;

        % Reset Qstay and Qleave values before running second block
        switch Model.Patch_sequence
            case 'Multiple'
                QStay = zeros(Task.Duration+1,Task.PatchNum);   % Q per patch type
            case 'Single'
                QStay = zeros(Task.Duration+1, 1);              % single Q over all patches
        end
        switch Model.Leave_sequence
            case 'Separate'
                QLeave = zeros(Task.TravelTime+1,1);            % all seconds in leave states are represented
            case 'Stateless'
                QLeave = 0;                                     % single leave state
        end

    end

end

%% Store variables
% Save QStay and QLeave
switch Model.Block_presentation
    case 'Combined'
        out.QStay = QStay(1:max(sum(QStay~=0)),:); % limit Q_stay tables - remove values not used
        out.QLeave = QLeave;
        out.NumObservations = NumObservations{BlockOrder(b)};
    case 'Separate'
        %Store the last block
        QStayBlock{BlockOrder(b)}= QStay(1:max(sum(QStay~=0)),:);
        QLeaveBlock{BlockOrder(b)}= QLeave;
        %Save
        out.QStay= QStayBlock;
        out.QLeave= QLeaveBlock;
        out.NumObservations = NumObservations;
end

if strcmpi(Model.Patch_sequence,'Multiple')
    out.NextPatch = NextPatch;
end

out.Q = Q;
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







