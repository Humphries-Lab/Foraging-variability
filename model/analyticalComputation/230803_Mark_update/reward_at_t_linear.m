function reward_in_patch_at_t = reward_at_t_linear(t,r0,alpha)

% REWARD_AT_T_linear linearly-decaying reward obtained in patch at time-step T
% R = REWARD_AT_T_LINEAR(T,RO,ALPHA) computes the reward R obtained in a patch
% at time-step T, given a linearly decaying amount of reward
% parameterised by initial reward in patch R0 and decay rate ALPHA
%
% R(T) = R0 - ALPHA*T
% Note: R(T) is set to 0 if above evaluates below 0
%
% Mark Humphries 28/7/2023

reward_in_patch_at_t = r0 - alpha*t;

reward_in_patch_at_t(reward_in_patch_at_t < 0) = 0;

