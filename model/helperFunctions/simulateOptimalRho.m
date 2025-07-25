function rho = simulateOptimalRho(g, rewardFunction, r0_list, patchOrder, travelTime, blockTime)
% simulateOptimalRR_timeLimited - Computes optimal average reward rate (rho) under a fixed time budget
%
% Inputs:
%   g               - decay rate
%   r0_list       - vector of initial reward rates for each patch type (e.g. [32.5, 45, 57.5])
%   patchOrder      - vector of patch type indices
%   travelTime      - travel time between patches
%   blockTime  - total duration of foraging block (in seconds)
%
% Output:
%   rho             - optimal average reward rate based on task parameters

% Initial guess and parameters
rho = 30;
rho_tol = 1e-4;
max_iter = 1000;

for iter = 1:max_iter
    totalReward = 0;
    cumulativeTime = 0;

    for i = 1:length(patchOrder)
        if cumulativeTime >= blockTime
            break;
        end

        patchType = patchOrder(i);
        r0 = r0_list(patchType);
        %patchTotalReward = 0;

        % Compute optimal leaving time
        if rho <= 0 || rho >= r0
            t_star = 0;
            reward = 0;
        else
            if strcmp(rewardFunction, 'exponential')
                t_star = log(rho / r0) / -g;
                reward = integral(@(x) r0 * exp(-g * x), 0, t_star);
            elseif strcmp(rewardFunction, 'linear')
                t_star = (rho - r0) / -g;
                reward = integral(@(x) max(0, r0 - g * x), 0, t_star); 

            end
        end

        patchTotalTime = t_star + travelTime;

        % If not enough time to complete this patch, truncate
        if cumulativeTime + patchTotalTime > blockTime
            remainingTime = blockTime - cumulativeTime;

            % Time in patch (minus travel)
            t_patch = max(remainingTime - travelTime, 0);
            %patchTotalReward = 0;
            if strcmp(rewardFunction, 'exponential')
                reward = integral(@(x) r0 * exp(-g * x), 0, t_patch);
            elseif strcmp(rewardFunction, 'linear')
                reward = integral(@(x) max(0, r0 - g * x), 0, t_patch);
            end

            cumulativeTime = blockTime;
            break;
        else
            totalReward = totalReward + reward;
            cumulativeTime = cumulativeTime + patchTotalTime;
        end
    end

    new_rho = totalReward / blockTime;

    if abs(new_rho - rho) < rho_tol
        break;
    end

    old_rho = rho;

    rho = new_rho;


end

end
