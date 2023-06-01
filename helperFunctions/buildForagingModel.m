function [sim_f, NLL_f, paramsIndex] = buildForagingModel(model, Env, subjData, startValues)

foragingModels = {
% Model N,       % simulation function                                % NLL% fitting function             % required parameters [alpha patch, alpha rho, beta]
    1, @(x)simulate_MVT_avgRR_RW_softmax(x,Env,startValues),   @(x0)NLL_MVT_avgRR_RW_softmax(x0, Env, subjData),   logical([0 1 1]);
    2, @(x)simulate_MVT_patchRR_RW_softmax(x,Env,startValues), @(x0)NLL_MVT_patchRR_RW_softmax(x0, Env, subjData), logical([1 0 1]);
    3, @(x)simulate_MVT_full_RW_softmax(x,Env,startValues),    @(x0)NLL_MVT_full_RW_softmax(x0, Env, subjData),    logical([1 1 1]);
    4, @(x)simulate_MVT_softmax(x,Env,startValues),            @(x0)NLL_MVT_softmax(x0, Env, subjData),            logical([0 0 1]);
    5, @(x)simulate_MVT_patchRR_RW_fixbeta(x,Env,startValues), @(x0)NLL_MVT_patchRR_RW_fixbeta(x0, Env, subjData), logical([1 0 0]);
};

% Assign the parameters index 

sim_f = foragingModels{model,2}; % simulation function 
NLL_f = foragingModels{model,3}; % NLL function 
paramsIndex = foragingModels{model, 4}; % which parameters the model needs