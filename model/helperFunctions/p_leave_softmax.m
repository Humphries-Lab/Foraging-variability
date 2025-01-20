function p = p_leave_softmax(reward_at_t,beta, bias, rho)

% P_LEAVE_SOFTMAX compute the probability of leaving the patch
% P = P_LEAVE_SOFTMAX(R,BETA) computes probability of leaving a patch given
% just-obtained reward R and softmax temperature BETA
%
% 31/7/2023: edit to handle R as a vector
% 08/1/2024: edit to add bias/intercept parameter (Emma S)
% 12/1/2024: edit to add rho for learning models (Emma S) 
% Mark H 5/5/2023

% softmax probability of leaving
p = (1 + exp(bias+beta .* (reward_at_t-rho))).^-1;
