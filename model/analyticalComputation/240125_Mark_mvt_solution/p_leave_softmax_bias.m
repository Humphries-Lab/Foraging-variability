function p = p_leave_softmax_bias(reward_at_t,beta,bias)

% P_LEAVE_SOFTMAX compute the probability of leaving the patch
% P = P_LEAVE_SOFTMAX(R,BETA,BIAS) computes probability of leaving a patch given
% just-obtained reward R and softmax temperature BETA and constant BIAS
%
% 31/7/2023: edit to handle R as a vector
% Mark H 24/1/2024

% softmax probability of leaving
p = (1 + exp(bias + beta .* reward_at_t)).^-1;
