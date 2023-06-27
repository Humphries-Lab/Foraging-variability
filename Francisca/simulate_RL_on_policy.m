function [output, Q_stay, Q_leave, rho, numObservations] = simulate_RL_on_policy(params, environ, block, model_table, model)


%% Set up

if block == 1
    patch_order = environ.patch_order.rich;
elseif block == 2
    patch_order = environ.patch_order.poor;
end

leave = 1;
stay = 2;

% reinforcement learning parameters
alpha_Q = params(1);
alpha_rho = params(2);
beta = params(3); % for SOFTMAX CHOICE. 
%epsilon = 0.1; % for EPSILON-GREEDY CHOICE. Can either fix or leave as
%free parameter 

leave_sequence = model_table.leave_sequence{model};
patch_sequence = model_table.patch_sequence{model};
next_patch_type = model_table.next_patch_type{model};

time_in_block = [1:environ.block_time+1]'; % plus 1 for the final Q update
Q = zeros(environ.block_time+1,2); % [Q(leave), Q(stay)]
PAction = zeros(environ.block_time+1,2); % probability of selecting leave or stay
Actions = zeros(environ.block_time+1,1);
patch_reward = zeros(environ.block_time+1,1); % discrete patch reward in each state
patchRR = zeros(environ.block_time+1,1); % patch reward rate
RPE = zeros(environ.block_time+1,1); % reward prediction error
rho = zeros(environ.block_time+1, 1); % estimated averageRR
PAction(1, :) = [0 1]; % set first probabilities (starting in patch)

numObservations = 0;

% initialise state values for Q(stay) based on model settings
switch patch_sequence
    case 'multiple'
        Q_stay = zeros(environ.block_time+1,3); % Q per patch type
    case 'single'
        Q_stay = zeros(environ.block_time+1, 1); % single Q over all patches
end
% initialise state values for Q(leave) based on model settings
switch leave_sequence
    case 'separate'
        Q_leave = zeros(environ.travel_time+1,1); % all seconds in leave states are represented
    case 'stateless'
        Q_leave = 0; % single leave state
end

% set the first actions for the first patch
PAction(1,:) = [0 1]; % 100% probability of staying
Actions(1) = 2; % set first action to stay

patch_number = 0; % start outside the patch
arrive = 1; % arriving into a new patch 
data = []; % initialise

% run model
for ii = 1:environ.block_time

    if Actions(ii) == stay % take action to stay

        if arrive % if arriving in the patch
            T = 1; % time in patch
            patch_number = patch_number + 1; % patch number increases
            patch_type = patch_order(patch_number); % select patch type
            if strcmp(patch_sequence, 'multiple')
                Q(ii,stay) = Q_stay(T, patch_type); % change value of first state based on what the patch actually is
            end
            arrive = 0;
        end

        % observe reward in CURRENT state
        patch_reward(ii) = environ.r(T,patch_type); % patch_reward depends on time in patch and patch type
        patchRR(ii) = patch_reward(ii)/environ.time_step; % patch reward

        % observe next state - log the next Q values to choose next action.
        switch patch_sequence
            case 'multiple'
                Q(ii+1,:) = [Q_leave(1), Q_stay(T+1,patch_type)]; % select Q_stay based on patch_type
            case 'single'
                Q(ii+1,:) = [Q_leave(1), Q_stay(T+1)]; % no patch type to select
        end

        % SOFTMAX CHOICE
        PAction(ii + 1,:) = exp(Q(ii+1, :) * beta) ./ sum(exp(beta*Q(ii+1, :)));
        Actions(ii+1) = discreteinvrnd(PAction(ii+1,:),1,1) ;    

        % EPSILON-GREEDY CHOICE 
%         if rand < epsilon % for some small probability epsilon
%             Actions(ii+1) = discreteinvrnd([0.5 0.5], 1, 1); % coin flip on action taken
%         else
%             Actions(ii+1) = max(Q(ii+1,:));
%         end
% one thing that could be helpful: when the agent first starts and Q-values
% are all 0, set the greedy action to 'stay'. This will help the model to
% learn initially, rather than get stuck leaving because it hasn't learned
% any patch values. 

        numObservations = numObservations+1;

        RPE(ii) = (patch_reward(ii)/57.5 - rho(ii)) + Q(ii+1, Actions(ii+1)) - Q(ii, Actions(ii));

        % update estimate of average RR
        rho(ii+1) = rho(ii) + alpha_rho * RPE(ii);

        % update Q table for staying
        switch patch_sequence
            case 'multiple'
                Q_stay(T,patch_type) = Q_stay(T,patch_type) + alpha_Q * RPE(ii);
            case 'single'
                Q_stay(T) = Q_stay(T) + alpha_Q * RPE(ii);
        end

        data = [data; patch_number, patch_type, T];

        T = T+1;

        if Actions(ii+1) == leave
            t = 1; % if next action is leave, then reset travel time counter
            % what does agent think the NEXT patch is?
            switch next_patch_type
                case 'previous'
                    next_patch = patch_type; % assumes next patch is the same as the current one
                case 'godlike'
                    next_patch = patch_order(patch_number+1);  % knows what the next patch is
                case 'stochastic'
                    if rand < 0.1 % for some % of the time, choose randomly
                        next_patch = randi([1 3]);
                    else % else assume most frequent patch
                        next_patch = mode(patch_order(1:patch_number));
                    end
            end
        end

    elseif Actions(ii) == leave % take action to leave
        T = nan; % not in a patch 
        switch patch_sequence
            case 'multiple'
                switch leave_sequence
                    case 'separate'
                        Q(ii+1,:) = [Q_leave(t+1), Q_stay(1,next_patch)]; % select Q_stay based on patch_type
                    case 'stateless'
                        Q(ii+1,:) = [Q_leave, Q_stay(1,next_patch)]; % select Q_stay based on patch_type
                end
            case 'single'
                switch leave_sequence
                    case 'separate'
                        Q(ii+1,:) = [Q_leave(t+1), Q_stay(1)];
                    case 'stateless'
                        Q(ii+1,:) = [Q_leave, Q_stay(1)];
                end
        end

        if t == environ.travel_time % if on the last second of travelling
            PAction(ii+1,:) = [0 1]; % force staying on next action
            Actions(ii+1) = stay;
            arrive = 1; % about to arrive in new patch

        elseif t < environ.travel_time % if still travelling
            PAction(ii+1,:) = [1 0]; % force leaving for duration of travel time
            Actions(ii+1) = leave;
        end

        % observe reward in CURRENT state
        patch_reward(ii) = 0; % not getting anything during travel
        patchRR(ii) = 0; % patch reward rate;

        RPE(ii) = (patch_reward(ii)/57.5 - rho(ii)) + Q(ii+1, Actions(ii+1)) - Q(ii, Actions(ii));

        % update estimate of average RR
        rho(ii+1) = rho(ii) + alpha_rho * RPE(ii);

        % update Q(leave) table
        switch leave_sequence
            case 'separate'
                Q_leave(t) = Q_leave(t) + alpha_Q * RPE(ii); % update Q-leave based on individual states
            case'stateless'
                Q_leave = Q_leave + alpha_Q * RPE(ii); % update all states of Q_leave
        end

        t = t+1;
        data = [data; patch_number, patch_type, T];

    end

end
data = [data; zeros([environ.block_time+1-ii,3])];
tmp = [PAction, Actions, data, time_in_block, Q, rho, patch_reward, patchRR, RPE];
rho = rho(1:ii);
output = array2table(tmp, 'VariableNames', {'p_leave', 'p_stay', 'a_selected', 'patch_number', 'patch_type', 'time_in_patch', 'time_in_block', 'Q_leave', 'Q_stay', 'rho', 'patch_reward', 'patchRR', 'RPE'});
end