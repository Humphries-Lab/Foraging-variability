function p = p_leave_softmax(reward_at_t,beta,epsilon)

% P_LEAVE_SOFTMAX compute the probability of leaving the patch
% P = P_LEAVE_SOFTMAX(R,BETA) computes probability of leaving a patch given
% just-obtained reward R and softmax temperature BETA
%
% Mark H 5/5/2023
% Emma S 16/05/23 - added 'lapse' component to convert to e-softmax 

% softmax probability of leaving
p = (1 + exp(beta * reward_at_t))^-1 * (1-2*epsilon) + epsilon;
