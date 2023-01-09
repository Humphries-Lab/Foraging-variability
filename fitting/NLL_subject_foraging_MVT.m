function [NegLogLikelihood,PAction, numObservations] = NLL_subject_foraging_MVT(params,Data,environ)

%% Set up

leave = 1;
stay = 2;

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
alpha_patchRR = params(1);
alpha_rho = params(2);
beta = params(3);

PAction = zeros(time_in_block,2); % probability of selecting leave or stay. Extra 1s for final update. 
estimated_patchRR = zeros(time_in_block,1); 
rho = zeros(time_in_block,1);

% initialise start values
PAction(1, :) = [0 1]; % set first probabilities (starting in patch)

arrive = 1;
LogLikelihood = 0; 
numObservations = 0; 

% run model
for ii = 1:time_in_block-1

    if Actions(ii) == stay % take action to stay

        if arrive % if arriving in the patch
            T = 1; % time in patch
            arrive = 0;
        end

        % observe reward and update estimates
            estimated_patchRR(ii+1) = estimated_patchRR(ii) + alpha_patchRR * (Rewards(ii)/57.5 - estimated_patchRR(ii));
            rho(ii+1) = rho(ii) + alpha_rho * (Rewards(ii)/57.5 - rho(ii));

            PAction(ii + 1,:) = exp([rho(ii+1), estimated_patchRR(ii+1)] * beta) ./ sum(exp(beta*[rho(ii+1), estimated_patchRR(ii+1)]));
            Pselected = PAction(ii+1, Actions(ii+1));
            numObservations = numObservations+1; % track how many data points used for Pselected 
            %(because it's not all actions due to forced leave/stay setup)
            %and because we remove data points past 600s. 

            if Pselected == 0  % prevent log(0) becoming infinity (and breaking fmincon)
                Pselected = eps(0); % make a very small number instead 
            end 

            LogLikelihood = LogLikelihood + log(Pselected);


        if Actions(ii+1) == leave
            t = 1; % if next action is leave, then reset travel time counter
        end

        T = T+1;

    elseif Actions(ii) == leave % take action to leave

        if t == environ.travel_time % if on the last second of travelling
            PAction(ii+1,:) = [0 1]; % force staying on next action
            arrive = 1; % about to arrive in new patch
        elseif t < environ.travel_time % if still travelling
            PAction(ii+1,:) = [1 0]; % force leaving for duration of travel time
        end

        % update estimate of average RR
        estimated_patchRR(ii+1) = estimated_patchRR(ii) + alpha_patchRR * (Rewards(ii)/57.5 - estimated_patchRR(ii));
        rho(ii+1) = rho(ii) + alpha_rho * (Rewards(ii)/57.5 - rho(ii));

        t = t+1;

    end
end

NegLogLikelihood = -LogLikelihood;  % return negative log likelihood: a positive number, which we seek to minimise

end
