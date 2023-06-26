function p = p_leave_epsilon(epsilon)

% P_LEAVE_EPSILON compute the probability of leaving the patch
% P = P_LEAVE_EPSILON(EPSILON) computes probability of leaving a patch given
% greedy parameter EPSILON
%
% Emma S 16/05/23

% probability of leaving

if rand < epsilon
    p = randi(0:1); % coin flip (stay/leave)
else
    p = 0; % assume stay as default action (because reward never < 0, so stay will always be max action) 
end
