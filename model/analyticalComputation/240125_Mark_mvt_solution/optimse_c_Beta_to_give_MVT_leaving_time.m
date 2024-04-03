% optimse (c,Beta) of bias model to find set of MVT leaving times
%
% Mark Humphries 25/1/24

clearvars; close all

% environment parameters for exponential decay of reward
r0 = [32.5, 45, 57.5];  % initial patch yield
R = 21.8678;  % average reward rate
alpha = 0.075;  % decay of reward

% solution parameters
n_steps = 1000;         % limit on time steps for computing expectation of leave time
params0 = [1,-10];    % initial guess values of parameters to intialise search (beta,bias)
lowerBounds = [0 -20];  % [beta, bias]
upperBounds = [1,0];  % arbitrary upper bound on beta to stop pathological behaviour

% don't allow bias: original beta-only model (can't fit shortest Expected
% leaving times)
%params0 = [1,0];    % initial guess values of parameters to intialise search (beta,bias)
%lowerBounds = [0 0];  % [beta, bias]
%upperBounds = [100,0];  % arbitrary upper bound on beta to stop pathological behaviour


% find solution...
f = @(x)compute_error(x,n_steps,r0,R,alpha);  % make anonymous function
[params_fitted,min_errorMVT] = fmincon(f,params0,[],[],[],[],lowerBounds,upperBounds);

% evaluate solution
for iPatch = 1:numel(r0)
    fittedExpectedTime(iPatch) = E_leave_bias(params_fitted(1),params_fitted(2),n_steps,r0(iPatch),alpha);
end

expectedMVT = -log(R./r0)*(1/alpha);

figure
plot(r0,expectedMVT,'.-'); hold on
plot(r0,fittedExpectedTime,'.-')
xlabel('Initial patch yield (r0)')
ylabel('Expected leaving time')
legend({'MVT','Model'},'Location','SouthEast')

%% functions
% error function to minimise 
function errorMVT = compute_error(x,n_steps,r0,R,alpha)
    beta = x(1);
    bias = x(2);

   expectedMVT = -log(R./r0)*(1/alpha);  % vector of eMVT errors

    for iPatch = 1:numel(r0)
        ExpectedTime(iPatch) = E_leave_bias(beta,bias,n_steps,r0(iPatch),alpha);
    end

    errorMVT = sum(abs(ExpectedTime - expectedMVT));  % try abs error first
end

% E(leave) function for bias model
function Expectation = E_leave_bias(beta,bias,n_steps,r0,alpha)
    for n = 1:n_steps
        r(n) = reward_at_t_exp(n,r0,alpha); % reward on that time-step
        p(n) = p_leave_softmax_bias(r(n),beta,bias);          % probability on that time-step
    
        % probability of leaving on that time-step is the product of
        % all probabilities of *not leaving* up till the current trial and the probability of
        % leaving on this trial
        if n == 1
            f_leave(n) = p(n);
        else
            f_leave(n) = prod(1-p(1:n-1)) * p(n); 
        end
    
    end
    Expectation = sum((1:n_steps).*f_leave);
end