% script to Monte Carlo generate leave-time distributions for different
% action selection policies - softmax and epsilon-softmax
% Sweep over averageRR, for fixed lambda value, for 3 different patch types (r0)
% 
% Mark Humphries 13 May 2023. 
% Emma Scholey updated - latest 12 July 2023

clearvars
close all;

number_samples = 250;  % number of Monte Carlo samples per parameter pair

rho = logspace(0,1.7,100);  % space of lambda values to calculate - based on range of participant fits
                        
t_max = 100;   % maximum time in patch (for explicit calculations)
t_step = 1;  % time-step at which to calculate estimates of E and VAR
alpha = 0.075;  % decay rate from LeHeron et al 202
r0 = [32.5 45 57.5];  % initial patch values from LeHeron et al 2020

n_steps = round(t_max / t_step);  % number of time-steps

lambda = 3;

%% MC sampling
t_series = (t_step* (1:n_steps))';  % the sequence of time-steps tested

sample_leave_time = zeros(numel(rho),numel(r0),number_samples);

for iL = 1:numel(rho)
    for iR = 1:numel(r0)
       %% calculate expected probability of leaving on each time-step
        f_leave_rich = zeros(n_steps,1); p_rich = zeros(n_steps,1); r = zeros(n_steps,1);
        f_leave_poor = zeros(n_steps,1); p_poor = zeros(n_steps,1);
        for n = 1:n_steps
            time_now = n * t_step;                          % what is actual time on this time-step?          
            r(n) = reward_at_t_exp(time_now,r0(iR),alpha);  % reward on that time-step
            beta_rich = lambda*1/(rho(iL)^0.2); % beta is some function of averageRR, based on rich environment
            beta_poor = lambda*1/(rho(iL)^0.2); % beta is some function of averageRR, based on poor environment

            p_rich(n) = p_leave_softmax(r(n),beta_rich);  % probability on that time-step
            p_poor(n) = p_leave_softmax(r(n),beta_poor);

            % probability of leaving on that time-step is the product of
            % all probabilities of *not leaving* up till the current trial and the probability of
            % leaving on this trial
            if n == 1
                f_leave_rich(n) = p_rich(n);
                f_leave_poor(n) = p_poor(n);

            else
                f_leave_rich(n) = prod(1-p_rich(1:n-1)) * p_rich(n); 
                f_leave_poor(n) = prod(1-p_poor(1:n-1)) * p_poor(n); 

            end
        end

        E_leave_rich(iL,iR) = sum(t_series.*f_leave_rich);
        E_leave_poor(iL,iR) = sum(t_series.*f_leave_poor);

        VAR_leave_rich(iL,iR) = sum((t_series - E_leave_rich(iL,iR)).^2 .* f_leave_rich); 
        VAR_leave_poor(iL,iR) = sum((t_series - E_leave_poor(iL,iR)).^2 .* f_leave_poor); 

    end
end

%% derived statistics from calculations
SD_leave_rich = sqrt(VAR_leave_rich);
SD_leave_poor = sqrt(VAR_leave_poor);

%% plot results
figure
semilogx(rho,E_leave_rich,'k'); hold on
xlabel('Rho: averageRR')
ylabel('Expected leaving time (s)')

figure
semilogx(rho,SD_leave_rich,'k'); hold on
xlabel('Rho: averageRR')
ylabel('SD of leaving time (s)')
