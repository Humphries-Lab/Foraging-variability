function [BIC, iBEST, BEST] = fitAll(Data, Env)

nStarts = 5; 
maxAlpha = [1 1];


% fit all of the models back to the simulated data of the model
[~, ~, BIC(1)] = fitM3_RLOn(Data,Env, nStarts, maxAlpha);
[~, ~, BIC(2)] = fitM4_RLOn(Data,Env, nStarts, maxAlpha);
[~, ~, BIC(3)] = fitM25_RLOn(Data,Env, nStarts, maxAlpha);
[~, ~, BIC(4)] = fitM26_RLOn(Data,Env, nStarts, maxAlpha);

% what's the best model for this data ?
[M, iBEST] = min(BIC);
BEST = BIC == M;
BEST = BEST / sum(BEST);