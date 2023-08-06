function p = p_leave_lapse(reward_at_t,beta, epsilon)

% P_LEAVE_LAPSE compute the probability of leaving the patch
% P = P_LEAVE_LAPSE(R,BETA, epsilon) computes probability of leaving a patch given
% just-obtained reward R, softmax temperature BETA, and lapse parameter
% EPSILON
%
% Emma S 16/05/23

% softmax probability of leaving
p = ((1 + exp(beta .* reward_at_t)).^-1) .* (1-2*epsilon) + epsilon;
