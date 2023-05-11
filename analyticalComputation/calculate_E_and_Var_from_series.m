% script to check calculation of E and Var of leaving time from series of
% expectations, where time t is a discrete variable
% 
% Note this is an approximation using discrete time-steps, like the
% simulated model. It thus predicts that beta also depends on the time-step
% of the simulation
% 
% Mark Humphries May 2023. 

clearvars
close all;

p_fixed = 0.1;

beta = logspace(-3,0);  % space of softmax temperatures to calculate; 
                        % maximum of beta=1 here, as actual temperature is
                        % beta * r0
                        
t_max = 100;   % maximum time in patch
t_step = 1;  % time-step at which to calculate estimates of E and VAR
alpha = 0.075;  % decay rate from LeHeron et al 202
r0 = [32.5 45 57.5];  % initial patch values from LeHeron et al 2020

n_steps = round(t_max / t_step);

%% validate against identical Bernoulli trials: 
% E and VAR should be from geometric distribution

% calculate expected probability of leaving on each time-step
f_leave = zeros(n_steps,1);
for n = 1:n_steps
    p = p_fixed;  % probability of event per time-step
    f_leave(n) = (1-p)^(n-1)*p;  % probability of leaving on that time-step
end

t_series = (1:n_steps)';
E_leave_fixed = sum(t_series.*f_leave)
geometric_dist_E = 1/p_fixed

VAR_leave_fixed = sum((t_series - E_leave_fixed).^2 .* f_leave)
geometric_dist_Var = (1-p_fixed) / p_fixed^2

%% test calculate for LeHeron task, with nonidentical Bernoulli trials
iB = 1;  % just use first beta
iR = 1;  % use first patch type
% calculate expected probability of leaving on each time-step
f_leave = zeros(n_steps,1); p = zeros(n_steps,1); r = zeros(n_steps,1);
for n = 1:n_steps
    time_now = n * t_step;                          % what is actual time on this time-step?
    r(n) = reward_at_t_exp(time_now,r0(iR),alpha); % reward on that time-step
    p(n) = p_leave_softmax(r(n),beta(iB));          % probability on that time-step
    f_leave(n) = n * (1-p(n))^(n-1)*p(n);         % probability of leaving on that time-step
end

t_series = (t_step* (1:n_steps))';  % the sequence of time-steps tested
E_leave = sum(t_series.*f_leave)
VAR_leave = sum((t_series - E_leave).^2 .* f_leave)

% does mean p approximate? No!
p_estimate = mean(p);
geometric_dist_E = 1/mean(p_estimate);
geometric_dist_Var = (1-p_estimate) / p_estimate^2;

%% now compute E and VAR across range of beta
t_series = (t_step* (1:n_steps))';  % the sequence of time-steps tested

for iB = 1:numel(beta)
    for iR = 1:numel(r0)
       % calculate expected probability of leaving on each time-step
        f_leave = zeros(n_steps,1); p = zeros(n_steps,1); r = zeros(n_steps,1);
        for n = 1:n_steps
            time_now = n * t_step;                          % what is actual time on this time-step?          
            r(n) = reward_at_t_exp(time_now,r0(iR),alpha);  % reward on that time-step
            p(n) = p_leave_softmax(r(n),beta(iB));          % probability on that time-step
            f_leave(n) = n * (1-p(n))^(n-1)*p(n);           % probability of leaving on that time-step
        end

        E_leave(iB,iR) = sum(t_series.*f_leave);
        VAR_leave(iB,iR) = sum((t_series - E_leave(iB,iR)).^2 .* f_leave); 
    end
end

% derived statistics
SD_leave = sqrt(VAR_leave);
CV_leave = E_leave ./ SD_leave;

%% plot results
figure
semilogx(beta,E_leave)
xlabel('Beta: higher = exploit')
ylabel('Expected leaving time (s)')

figure
semilogx(beta,VAR_leave)
xlabel('Beta: higher = exploit')
ylabel('Variance of leaving time (s)')

figure
semilogx(beta,SD_leave)
xlabel('Beta: higher = exploit')
ylabel('SD of leaving time (s)')

figure
semilogx(beta,CV_leave)
xlabel('Beta: higher = exploit')
ylabel('Coeff Var of leaving time (s)')

% plot examples like LeHeron task performance
% roughly from Fig 2A:
% rich: [11 16 21]
% poor: [15 20 25]

% find index that agrees with middle-valued patch (as it occured at same
% rate in both environments)
ixRich = find(E_leave(:,2) > 16,1,'first');  
ixPoor = find(E_leave(:,2) > 20,1,'first'); 

figure 
plot(E_leave(ixRich,:),'.-','Color',[0.7 0.3 0.3]); hold on
plot(E_leave(ixPoor,:),'.-','Color',[0.3 0.3 0.7]);
set(gca,'XLim',[0 4],'XTick',[ 1 2 3],'XTickLabel',{'Low','Medium','High'})
ylabel('Expected leaving time (s)')

beta_rich = beta(ixRich)
beta_poor = beta(ixPoor)

%% question: which is bigger - change in SD of patch-type between environments, 
% or between patch types in same environment?

% do for LeHeron estimates here - this is a more general question 

change_in_SD_between_environments = SD_leave(ixPoor,:) - SD_leave(ixRich,:)
change_in_SD_within_rich_environment = SD_leave(ixRich,:) - SD_leave(ixRich,:)'
change_in_SD_within_poor_environment = SD_leave(ixPoor,:) - SD_leave(ixPoor,:)'
% ANSWER: SDs between low and high should differ more within environment
% than patch SDs do between environment


