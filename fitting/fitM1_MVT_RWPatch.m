function [minNLL, minNLLFitParams, BIC, AIC, FitParams, NLLEval, StartParams] = fitM1_MVT_RWPatch(SubjData, Env, nStarts, maxAlpha)

% fmincon options
%lowerBounds = [0,0];  % [alpha_Q, alpha_rho, beta] % parameter bounds
%upperBounds = [1,100];  % arbitrary upper bound on beta to stop pathological behaviour
lowerBounds = [0];  % [alpha_Q, alpha_rho, beta] % parameter bounds
upperBounds = [1];  % arbitrary upper bound on beta to stop pathological behaviour
options = optimoptions('fmincon','Display','none'); % don't display
numParams = length(lowerBounds); 

FitParams = zeros([nStarts, numParams]);
StartParams = zeros([nStarts, numParams]);
NLLEval = zeros([nStarts,1]);

%% Run fitting
f = @(x)NLL_M1_MVT_RWPatch(x, Env, SubjData);

parfor ii = 1:nStarts
    params0 = [rand]; % choose start parameters
    StartParams(ii, :) = params0;
    [FitParams(ii,:),NLLEval(ii)] = fmincon(f,params0,[],[],[],[],lowerBounds,upperBounds, [], options);
end

%% Find the best fitting parameter values 
minNLL = min(NLLEval);   % minimum negative log likelihood over all starting positions
ix = find(minNLL == NLLEval);    % indices of location of minimum, to find the corresponding best fit parameters

% get the corresponding fitted parameter values. Note: if more than one starting location converges on same parameters then will have same NLL, so use first
minNLLFitParams = FitParams(ix(1),:);

% get the corresponding hessian
% minNLLFitParamsHessian = hessian(:,:,ix);

% calculate standard error on the parameter fits using the hessian matrix
% The diagonal terms of H−1 correspond to variances for each parameter separately, and their square roots measure one standard error on the parameter. - Daw 2011
% tmp = sqrt(diag(inv(minNLLFitParamsHessian)))';
% tmp(~isreal(tmp)) = nan;
% minNLLFitParamsSE = tmp;

%% Calculate BIC/AIC
sum(SubjData == 2); % get the number of data points that were fit to (equivalent to how many stay states there are e.g. states in the patch where they decide whether to stay or leave)

% BIC and AIC
BIC = numParams * log(numObv) + 2*minNLL;
AIC = 2/numObv * minNLL + 2 * numParams/numObv;
end

