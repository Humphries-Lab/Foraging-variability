function reward_in_patch_at_t = reward_at_t_exp(t,r0,alpha)

% REWARD_AT_T_EXP exponentially-decaying reward obtained in patch at time-step T
% R = REWARD_AT_T_EXP(T,RO,ALPHA) computes the reward R obtained in a patch
% at time-step T, given an exponentially decaying amount of reward
% parameterised by initial reward in patch R0 and decay rate ALPHA
%
% Mark Humphries 9/5/2023

reward_in_patch_at_t = r0*exp(-alpha*t);
