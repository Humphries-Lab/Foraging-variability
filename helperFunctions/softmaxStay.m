function p = softmaxStay(beta, valueStay, valueLeave)

% softmaxStay computes the probability of staying in the current patch
% P = softmaxStay(BETA, VALUE_STAY, VALUE_LEAVE) computes probability of
% staying in patch given the value of staying (either Q value or patch
% reward rate), value of leaving (either Q value or average reward rate)
%  and softmax temperature BETA
%
% Emma S 15/05/23

p = (1 + exp(- beta * (valueStay-valueLeave)))^-1;
