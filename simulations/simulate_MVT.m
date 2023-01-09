function [output, numObservations] = simulate_MVT(params, environ, block, model_table, model)

%% Set up

if block == 1
    patch_order = environ.patch_order.rich;
    averageRR = environ.averageRR.rich;
elseif block == 2
    patch_order = environ.patch_order.poor;
    averageRR = environ.averageRR.poor;
end

leave = 1;
stay = 2;

% reinforcement learning parameters
alpha_patchRR = params(1);
alpha_rho = params(2);
beta = params(3);

time_in_block = [1:environ.block_time+1]'; % plus 1 for the final Q update
PAction = zeros(environ.block_time+1,2); % probability of selecting leave or stay
Actions = zeros(environ.block_time+1,1);
patch_reward = zeros(environ.block_time+1,1); % reward prediction error
patchRR = zeros(environ.block_time+1,1); % real patch reward rate
estimated_patchRR = zeros(environ.block_time+1,1); % estimated patch reward rate
rho = zeros(environ.block_time+1,1); % estimated average reward rate

PAction(1, :) = [0 1]; % set first probabilities (starting in patch)
Actions(1) = 2; % set first action to stay

patch_number = 0;
arrive = 1;
data = [];
numObservations = 0; 

% run model
for ii = 1:environ.block_time

    if Actions(ii) == stay % take action to stay

        if arrive % if arriving in the patch
            T = 1; % time in patch
            patch_number = patch_number + 1; % patch number increases
            patch_type = patch_order(patch_number); % select patch type
            arrive = 0;
        end

        % observe reward in CURRENT state
        patch_reward(ii) = environ.r(T,patch_type); % patch_reward depends on time in patch and patch type
        patchRR(ii) = patch_reward(ii)/environ.time_step; % patch reward

        if strcmp(model_table.policy{model}, 'MVT_learning')
            % update agent estimates
            estimated_patchRR(ii+1) = estimated_patchRR(ii) + alpha_patchRR * (patch_reward(ii)/57.5 - estimated_patchRR(ii));
            rho(ii+1) = rho(ii) + alpha_rho * (patch_reward(ii)/57.5 - rho(ii));

            PAction(ii + 1,:) = exp([rho(ii+1), estimated_patchRR(ii+1)] * beta) ./ sum(exp(beta*[rho(ii+1), estimated_patchRR(ii+1)]));
            Actions(ii+1) = (rand(1) < PAction(ii + 1,2)) + 1;     % make a choice of next action: choose stim 1 if random number is in the [0 p(1)] interval, and 2 otherwise
            numObservations = numObservations+1; % track how many data points used for Pselected 

        elseif strcmp(model_table.policy{model}, 'MVT')
            Actions(ii+1) = (patchRR(ii) > averageRR) + 1;
        end

        data = [data; patch_number, patch_type, T];

        if Actions(ii+1) == leave
            t = 1; % if next action is leave, then reset travel time counter
        end

        T = T+1;

    elseif Actions(ii) == leave % take action to leave

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

        % update estimate of average RR
        estimated_patchRR(ii+1) = estimated_patchRR(ii) + alpha_patchRR * (patch_reward(ii)/57.5 - estimated_patchRR(ii));
        rho(ii+1) = rho(ii) + alpha_rho * (patch_reward(ii)/57.5 - rho(ii));

        t = t+1;
        T = nan; 
        data = [data; patch_number, patch_type, T];

    end
end

data = [data; NaN NaN NaN]; % add on final row for updatd parameters in ii+1
output = [PAction, Actions, time_in_block, data, estimated_patchRR, rho, patch_reward, patchRR];
output = array2table(output, 'VariableNames', {'p_leave', 'p_stay', 'a_selected', 'time_in_block', 'patch_number', 'patch_type', 'time_in_patch', 'estimated_patchRR', 'rho', 'patch_reward', 'patchRR'});
end
