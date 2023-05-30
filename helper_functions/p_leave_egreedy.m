function p = p_leave_egreedy(reward_at_t, epsilon)

p_leave = 0.2; % 2 is number of actions 
PAction(ii+1, ExploreAction) = (epsilon/2);

if rand < epsilon % for some small probability epsilon
    Action(ii+1) = discreteinvrnd([0.5 0.5], 1, 1); % coin flip on action taken
else
    Action(ii+1) = GreedyAction;
end