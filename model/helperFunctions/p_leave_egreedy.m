function p = p_leave_egreedy(reward_at_t,avgRR_at_t,epsilon)

% P_LEAVE_EPSILON compute the probability of leaving the patch
% P = P_LEAVE_EPSILON(EPSILON) computes probability of leaving a patch given
% greedy parameter EPSILON given just-obtained reward and avgRR

% Emma S 16/05/23

% probability of leaving

p = [0 0];
if rand < epsilon
    p(1) = randi(0:1); % coin flip (stay/leave)
    p(2) = 1- p(1); 
else
    max_action = [avgRR_at_t,reward_at_t] == max(avgRR_at_t,reward_at_t);
    p = double(max_action); % assume stay as default action
end