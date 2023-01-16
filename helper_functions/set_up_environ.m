%% Set up foraging environment according to Le Heron et al. (2020) 

% task parameters
environ.travel_time = 6; % delay between patches
environ.time_step = 1;  % T, seconds to update time_in_block and time_in_patch
environ.block_time = 1200;  % changing this affects fitting, so remember to change it back if playing around with simulations
patch_yields = [32.5, 45, 57.5]; % S in LeHeron paper

for j=1:length(patch_yields)
    fun = @(T) patch_yields(j)*exp(-T*0.075); %Exponential function currently
    rr(:,j) = fun(0:50);
    k = 0:environ.time_step:environ.block_time; %this is set for a maximum of block_time secs in the patch, with loop time according to time_step
    for i = 1:length(k)
        r(i,j) = integral(fun,k(i),k(i)+environ.time_step); % This creates matrix 'q' with each cell reward earned in time_step epoch at that reward rate
    end
end
clear i j k
% specify which patch the agent enters. This trial order was set by
% LeHeron for cabergoline data, where every 10 patches leads to correct
% proportions of low, medium or high patches based on the environment type
% (but randomised for PD data, where proportions are probabilities of encountering different patch types)

environ.patch_order.rich = repmat([3,2,3,1,3,2,1,3,3,2,1,3,2,1,3,3,2,3,3,2,3,3,1,3,3,2,3,1,2,2,3,2,2,1,3,3,2,3,1,3,3,3,2,2,3,3,3,1,1,2,3,2,2,3,3,1,1,3,2,3,3,2,3,1,1,3,2,2,3,3,3,3,3,1,3,2,3,2,2,1,1,3,2,3,3,1,2,2,3,3,2,3,1,3,3,3,1,2,2,3], [1,5]);
environ.patch_order.poor = repmat([1,1,2,1,3,2,2,1,3,1,1,2,3,1,2,1,1,2,3,1,1,1,2,3,1,2,3,2,1,1,3,1,1,2,2,1,2,3,1,1,1,2,3,1,1,2,1,3,1,2,3,2,2,1,2,1,1,1,3,1,1,1,2,1,2,1,1,2,3,3,1,3,1,3,2,2,2,1,1,1,2,3,1,3,1,2,1,2,1,1,2,1,2,3,2,3,1,1,1,1], [1,5]);

environ.averageRR.rich = 21.8678;
environ.averageRR.poor = 18.5632;

environ.r = r;