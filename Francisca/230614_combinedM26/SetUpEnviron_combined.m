%% Set up foraging environment according to Le Heron et al. (2020) 

% task parameters
Env.TravelTime = 6; % delay between patches
Env.TimeStep = 1;  % T, seconds to update time_in_block and time_in_patch
Env.BlockTime = 600;  % changing this affects fitting, so remember to change it back if playing around with simulations
Env.TaskTime = Env.BlockTime*2; % total time in task (2 blocks)
Env.InitialYield = [32.5, 45, 57.5]; % S in LeHeron paper

for j=1:length(Env.InitialYield) % for each yield type
    fun = @(T) Env.InitialYield(j)*exp(-T*0.075); %Exponential function
    Env.RR(:,j) = fun(0:Env.TaskTime); % for theoretical maximum of Task Time in the patch
    k = 0:Env.TimeStep:Env.TaskTime; %this is set for a maximum of block_time secs in the patch, with loop time according to time_step
    for i = 1:length(k)
        Env.R(i,j) = integral(fun,k(i),k(i)+Env.TimeStep); % This creates matrix with each cell reward earned in time_step epoch at that reward rate
    end
end
clear i j k
% specify which patch the agent enters. This trial order was set by
% LeHeron for cabergoline data, where every 10 patches leads to correct
% proportions of low, medium or high patches based on the environment type
% (but randomised for PD data, where proportions are probabilities of encountering different patch types)

% patch order in rich environment
Env.PatchOrder{1} = repmat([3,2,3,1,3,2,1,3,3,2,1,3,2,1,3,3,2,3,3,2,3,3,1,3,3,2,3,1,2,2,3,2,2,1,3,3,2,3,1,3,3,3,2,2,3,3,3,1,1,2,3,2,2,3,3,1,1,3,2,3,3,2,3,1,1,3,2,2,3,3,3,3,3,1,3,2,3,2,2,1,1,3,2,3,3,1,2,2,3,3,2,3,1,3,3,3,1,2,2,3], [1,5]);
% patch order in poor environment
Env.PatchOrder{2} = repmat([1,1,2,1,3,2,2,1,3,1,1,2,3,1,2,1,1,2,3,1,1,1,2,3,1,2,3,2,1,1,3,1,1,2,2,1,2,3,1,1,1,2,3,1,1,2,1,3,1,2,3,2,2,1,2,1,1,1,3,1,1,1,2,1,2,1,1,2,3,3,1,3,1,3,2,2,2,1,1,1,2,3,1,3,1,2,1,2,1,1,2,1,2,3,2,3,1,1,1,1], [1,5]);

Env.OptimalAverageRR{1} = 21.8678; % optimal average RR in rich environment
Env.OptimalAverageRR{2} = 18.5632; % optimal average RR in poor environment

