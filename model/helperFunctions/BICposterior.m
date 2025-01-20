function p = BICposterior(BICs)

% BICPOSTERIOR compute approximate posterior model probability from BIC
% P = BICPOSTERIOR(B) given the BIC scores of models 1, 2, ..., k in the
% k-length array B, returns P the approximate posterior probability of each
% model
%
% Wagenmakers, E.-J. (2007) A practical solution to the pervasive problems of p values. 
% Psychon Bull Rev, 14, 779-804 
% Eq 11 of that paper
% 
% Mark Humphries 6/6/23

expBICs = exp(vpa(-0.5.*BICs));

p = expBICs ./ sum(expBICs,2);

p = double(p);