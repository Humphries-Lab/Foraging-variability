function [NegLogLikelihood,rho,Q_stay,Q_leave,PAction,numObservations] = NLL_subject_foraging_RL_off(params,Data,patch_order,environ,model_table,model)

%% Set up
leave = 1;
stay = 2;

% get model settings
leave_sequence = model_table.leave_sequence{model};
patch_sequence = model_table.patch_sequence{model};
next_patch_type = model_table.next_patch_type{model};
off_policy_options = model_table.off_policy_options{model};

% how long was their foraging task?
if size(Data(:,1),1) <= 600 % if under 600s (e.g. excluded some trials)
    time_in_block = size(Data(:,1),1); % set time_in_block to their actual time
    Actions = Data(:,2);
    Rewards = Data(:,3);
elseif size(Data(:,1),1) > 600 % if over 600s (e.g. they stayed in the last patch past 600s)
    time_in_block = environ.block_time; % limit to max block_time (600s);
    Actions = Data(1:environ.block_time,2); % bound their data to max block_time
    Rewards = Data(1:environ.block_time,3);
end

% reinforcement learning parameters
alpha_Q = params(1);
alpha_rho = params(2);
beta = params(3);

PAction = zeros(time_in_block,2); % probability of selecting leave or stay. Extra 1s for final update.
rho = zeros(time_in_block,1);
RPE = zeros(time_in_block,1); % reward prediction error
Q = zeros(time_in_block,2); % [Q(leave), Q(stay)]

% initialise state values for Q(stay) based on model settings
if strcmp(patch_sequence, 'multiple')
    Q_stay = zeros(environ.block_time,3); % Q per patch type
elseif strcmp(patch_sequence, 'single')
    Q_stay = zeros(environ.block_time, 1);
end

% initialise state values for Q(leave) based on model settings
if strcmp(leave_sequence, 'separate')
    Q_leave = zeros(environ.travel_time+1,1);
elseif strcmp(leave_sequence, 'stateless')
    Q_leave = 0;
end

PAction(1, :) = [0 1]; % set first probabilities (starting in patch)

LogLikelihood = 0;
numObservations = 0;

patch_number = 0;
arrive = 1;

% run model
for ii = 1:time_in_block-1

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

        % observe next state - log the next Q values to choose next action.
        if strcmp(patch_sequence, 'multiple')
            Q(ii+1,:) = [Q_leave(1), Q_stay(T+1,patch_type)]; % select Q_stay based on patch_type
        elseif strcmp(patch_sequence, 'single')
            Q(ii+1,:) = [Q_leave(1), Q_stay(T+1)]; % no patch type to select
        end

        RPE(ii) = (Rewards(ii)/57.5 - rho(ii)) + max(Q(ii+1, :)) - Q(ii, Actions(ii));

        % update estimate of average RR
        rho(ii+1) = rho(ii) + alpha_rho * RPE(ii);

        % update Q table for staying
        if strcmp(patch_sequence, 'multiple')
            Q_stay(T,patch_type) = Q_stay(T,patch_type) + alpha_Q * RPE(ii);
        elseif strcmp(patch_sequence, 'single')
            Q_stay(T) = Q_stay(T) + alpha_Q * RPE(ii);
        end

        softmax = exp(Q(ii+1, :) * beta); % softmax function
        softmax(softmax==inf) = realmax/2; % replace any infs with just a large number
        softmax(softmax==0) = eps(0); % replace any 0's with just a tiny number to prevent NaN emerging in softmax calculation
        PAction(ii + 1,:) = softmax ./ sum(softmax);

        Pselected = PAction(ii+1, Actions(ii+1));
        numObservations = numObservations+1;

        if Pselected == 0  % prevent log(0) becoming infinity (and breaking fmincon)
            Pselected = eps(0); % make a very small number instead
        end

        LogLikelihood = LogLikelihood + log(Pselected);

        T = T+1;

        if Actions(ii+1) == leave
            t = 1; % if next action is leave, then reset travel time counter
            T_at_leaving = T; % log state at which they left

            % what does agent think the NEXT patch is?
            if strcmp(next_patch_type, 'previous')
                next_patch = patch_type; % assumes next patch is the same as the current one
            elseif strcmp(next_patch_type, 'godlike')
                next_patch = patch_order(patch_number+1);  % knows what the next patch is
            elseif strcmp(next_patch_type, 'stochastic')
                if rand < 0.1 % for some % of the time, choose randomly
                    next_patch = randi([1 3]);
                else % else assume most frequent patch
                    next_patch = mode(patch_order(1:patch_number));
                end
            end
        end

    elseif Actions(ii) == leave % take action to leave

        T = nan; 

 switch patch_sequence
            case 'multiple'
                switch leave_sequence
                    case 'separate'
                        switch off_policy_options
                            case 'anticipate'
                                Q(ii+1,:) = [Q_leave(t+1), Q_stay(1,next_patch)]; % select Q_stay based on patch_type
                            case 'regret'
                                Q(ii+1,:) = [Q_leave(t+1), Q_stay(T_at_leaving,next_patch)];
                            case  'dynamic_regret'
                                Q(ii+1,:) = [Q_leave(t+1), Q_stay(T_at_leaving+t,next_patch)];
                        end
                    case 'stateless'
                        switch off_policy_options
                            case 'anticipate'
                                Q(ii+1,:) = [Q_leave, Q_stay(1,next_patch)];
                            case 'regret'
                                Q(ii+1,:) = [Q_leave, Q_stay(T_at_leaving,next_patch)];
                            case 'dynamic_regret'
                                Q(ii+1,:) = [Q_leave, Q_stay(T_at_leaving+t,next_patch)];
                        end
                end

            case 'single'
                switch leave_sequence
                    case 'separate'
                        switch off_policy_options
                            case 'anticipate'
                                Q(ii+1,:) = [Q_leave(t+1), Q_stay(1)]; % select Q_stay based on patch_type
                            case 'regret'
                                Q(ii+1,:) = [Q_leave(t+1), Q_stay(T_at_leaving)];
                            case 'dynamic_regret'
                                Q(ii+1,:) = [Q_leave(t+1), Q_stay(T_at_leaving+t)];
                        end

                    case 'stateless'
                        switch off_policy_options
                            case 'anticipate'
                                Q(ii+1,:) = [Q_leave, Q_stay(1)];
                            case 'regret'
                                Q(ii+1,:) = [Q_leave, Q_stay(T_at_leaving)];
                            case 'dynamic_regret'
                                Q(ii+1,:) = [Q_leave, Q_stay(T_at_leaving+t)];
                        end
                end
        end

        if t == environ.travel_time % if on the last second of travelling
            PAction(ii+1,:) = [0 1]; % force staying on next action
            arrive = 1; % about to arrive in new patch

        elseif t < environ.travel_time % if still travelling
            PAction(ii+1,:) = [1 0]; % force leaving for duration of travel time
        end

        RPE(ii) = (Rewards(ii)/57.5 - rho(ii)) + max(Q(ii+1, :)) - Q(ii, Actions(ii));

        % update estimate of average RR
        rho(ii+1) = rho(ii) + alpha_rho * RPE(ii);

        % update Q(leave) table
        if strcmp(leave_sequence, 'separate')
            Q_leave(t) = Q_leave(t) + alpha_Q * RPE(ii); % update Q-leave based on individual states
        elseif strcmp(leave_sequence, 'stateless')
            Q_leave = Q_leave + alpha_Q * RPE(ii); % update all states of Q_leave
        end

        t = t+1;
    end
end

NegLogLikelihood = -LogLikelihood;  % return negative log likelihood: a positive number, which we seek to minimise

end
