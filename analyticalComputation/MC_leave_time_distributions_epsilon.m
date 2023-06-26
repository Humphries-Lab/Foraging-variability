%MC_leave_time_distributions_epsilon

% script to Monte Carlo generate leave-time distriEutions for lapse policy
% Sweep over epsilon parameter, for 3 different patch types (r0). 
% Can look at two different choice rules - e-greedy and e-softmax. 
% 
% Emma Scholey 16 May 2023. 
% adapted from MC_leave_time_distributions from MH

clearvars
close all;

number_samples = 1000;  % number of Monte Carlo samples per parameter pair

beta = 0.25; % fix at 'best fit' for now - for epsilon-softmax model

epsilon = logspace(-3,0); % for lapse model
                        
t_max = 100;   % maximum time in patch (for explicit calculations)
t_step = 1;  % time-step at which to calculate estimates of E and VAR
alpha = 0.075;  % decay rate from LeHeron et al 202
r0 = [32.5 45 57.5];  % initial patch values from LeHeron et al 2020

n_steps = round(t_max / t_step);  % number of time-steps


%% MC sampling
t_series = (t_step* (1:n_steps))';  % the sequence of time-steps tested

sample_leave_time = zeros(numel(epsilon),numel(r0),number_samples);
E_leave = zeros(numel(epsilon),numel(r0)); VAR_leave = zeros(numel(epsilon),numel(r0));

for iE = 1:numel(epsilon)
    for iR = 1:numel(r0)
       %% calculate expected probability of leaving on each time-step
        f_leave = zeros(n_steps,1); p = zeros(n_steps,1); r = zeros(n_steps,1);
        for n = 1:n_steps
            time_now = n * t_step;                          % what is actual time on this time-step?          
            r(n) = reward_at_t_exp(time_now,r0(iR),alpha);  % reward on that time-step
            %p(n) = p_leave_lapse(r(n),beta, epsilon(iE));          % probability on that time-step
            p(n) = p_leave_epsilon(epsilon(iE));                % probability on that time-step - no function of reward

            % THIS IS NOT RIGHT: IT IS NOT 1-p(current step)^n-1: it is
            % cumulative 1-p of all previous failures to leave!
            % f_leave(n) = n * (1-p(n))^(n-1)*p(n);           % probability of leaving on that time-step
            
            % probability of leaving on that time-step is the product of
            % all probabilities of *not leaving* up till the current trial and the probability of
            % leaving on this trial
            if n == 1
                f_leave(n) = p(n);
            else
                f_leave(n) = prod(1-p(1:n-1)) * p(n); 
            end
        end

        E_leave(iE,iR) = sum(t_series.*f_leave);
        VAR_leave(iE,iR) = sum((t_series - E_leave(iE,iR)).^2 .* f_leave); 


        %% Monte Carlo sampling of distriEution
        for iSample = 1:number_samples
            blnLeave = 0; n = 1;
            while ~blnLeave
                time_now = n * t_step;                          % what is actual time on this time-step?          
                r = reward_at_t_exp(time_now,r0(iR),alpha);  % reward on that time-step
                %p = p_leave_lapse(r,beta, epsilon(iE));          % probability on that time-step
                 p = p_leave_epsilon(epsilon(iE));                % probability on that time-step - no function of reward

                blnLeave = rand <= p;
                n = n + 1;
            end
            sample_leave_time(iE,iR,iSample) = time_now;
        end
    end
end

%% derived statistics from calculations
SD_leave = sqrt(VAR_leave);
CV_leave = E_leave ./ SD_leave;

%% statistics from MC simulations
E_MC_leave = mean(sample_leave_time,3);
SD_MC_leave = std(sample_leave_time,[],3); 

%% plot results
figure
semilogx(epsilon,E_leave,'k'); hold on
semilogx(epsilon,E_MC_leave); hold on
xlabel('epsilon: higher = more noisy')
ylabel('Expected leaving time (s)')

figure
semilogx(epsilon,SD_leave,'k'); hold on
semilogx(epsilon,SD_MC_leave);
xlabel('epsilon: higher = more noisy')
ylabel('SD of leaving time (s)')

figure
semilogx(epsilon,CV_leave)
xlabel('epsilon: higher = more noisy')
ylabel('Coeff Var of leaving time (s)')

% plot examples like LeHeron task performance
% roughly from Fig 2A:
% rich: [11 16 21]
% poor: [15 20 25]

% find index that agrees with middle-valued patch (as it occured at same
% rate in both environments)
%ixRich = find(E_leave(:,2) > 16,1,'last');  
%ixPoor = find(E_leave(:,2) > 20,1,'last'); 
ixRich = find(E_MC_leave(:,2) > 16,1,'last');  % use MC version for e-greedy 
ixPoor = find(E_MC_leave(:,2) > 20,1,'last'); % use MC version for e-greedy
% 
figure 
plot(E_MC_leave(ixRich,:),'.-','Color',[0.7 0.3 0.3], 'LineWidth', 2); hold on
plot(E_MC_leave(ixPoor,:),'.-','Color',[0.3 0.3 0.7], 'LineWidth', 2);
set(gca,'XLim',[0 4],'XTick',[ 1 2 3],'XTickLabel',{'Low','Medium','High'})
ylabel('Expected leaving time (s)')

epsilon_rich = epsilon(ixRich)
epsilon_poor = epsilon(ixPoor)


%% question: which is bigger - change in SD of patch-type between environments, 
% or between patch types in same environment?

% do for LeHeron estimates here - this is a more general question 

change_in_SD_between_environments = SD_MC_leave(ixPoor,:) - SD_MC_leave(ixRich,:)
change_in_SD_within_rich_environment = SD_MC_leave(ixRich,:) - SD_MC_leave(ixRich,:)'
change_in_SD_within_poor_environment = SD_MC_leave(ixPoor,:) - SD_MC_leave(ixPoor,:)'
% ANSWER: SDs between low and high should differ more within environment
% than patch SDs do between environment

% E and SD at participant fit 
E_MC_leave(ixPoor,:)
E_MC_leave(ixRich,:)

SD_leave_mean(1,:) = SD_MC_leave(ixRich,:);
SD_leave_mean(2,:) = SD_MC_leave(ixPoor,:);

 %% compare to participant data

% initialise
numBlocks = 2;
numPatches = 3;

load t_young.mat % load summarised subject data 
numSubjects = numel(unique(t_young.subj));

% test analytical predictions 

SD_LT = zeros([numSubjects, numPatches, numBlocks]); % each row is subject SD, each column patch type, 3d is block type

for iSubj = 1:numSubjects % for each subject
    for iPatch = 1:numPatches % each patch
        for iBlock = 1:numBlocks % each block
            SD_LT(iSubj,iPatch,iBlock) = std(t_young.leaveT(t_young.subj == iSubj & t_young.env == iBlock & t_young.patch == iPatch));
        end
    end
end

% between environments 

real_change_in_SD_between_environments = mean(SD_LT(:,:,2) - SD_LT(:,:,1)); % poor - rich

% within environments 

real_change_in_SD_within_rich = mean(SD_LT(:,:,1)) - mean(SD_LT(:,:,1))';
real_change_in_SD_within_poor = mean(SD_LT(:,:,2)) - mean(SD_LT(:,:,2))';

%mean SDs
SD_LT_mean = squeeze(mean(SD_LT))';

% plot model vs participant SDs across patches and environment
figure

plot(SD_MC_leave(ixRich,:),'.--','Color',[0.7 0.3 0.3], 'LineWidth', 2, 'MarkerSize', 15); hold on
plot(SD_MC_leave(ixPoor,:),'.--','Color',[0.3 0.3 0.7], 'LineWidth', 2,'MarkerSize', 15);
plot(SD_LT_mean(1,:),'.-','Color',[0.7 0.3 0.3], 'LineWidth', 2,'MarkerSize', 15) % rich
plot(SD_LT_mean(2,:),'.-','Color',[0.3 0.3 0.7], 'LineWidth', 2,'MarkerSize', 15) % poor
scatter([1 2 3], SD_LT(:,:,1),[],[0.7 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.1) % rich
scatter([1 2 3], SD_LT(:,:,2),[],[0.3 0.3 0.7], 'filled', 'MarkerFaceAlpha', 0.1) % poor

set(gca,'XLim',[0 4],'XTick',[ 1 2 3],'XTickLabel',{'Low','Medium','High'})
%set(gca,'YLim', [0 12])
ylabel('SD of leaving time (s)')

% plot between environment differences 
figure

plot(change_in_SD_between_environments,'.--', 'Color', 'black', 'LineWidth', 2, 'MarkerSize', 15); hold on
plot(real_change_in_SD_between_environments,'.-', 'Color', 'black', 'LineWidth', 2,'MarkerSize', 15) % rich
scatter([1 2 3], SD_LT(:,:,2) - SD_LT(:,:,1), 'black',  'filled', 'MarkerFaceAlpha', 0.1) % rich

set(gca,'XLim',[0 4],'XTick',[ 1 2 3],'XTickLabel',{'Low','Medium','High'})
ylabel('Change in SD leaving time (s) (poor - rich)')


% plot within environment differences 
figure
change_rich = nonzeros(triu(change_in_SD_within_rich_environment))';
change_poor = nonzeros(triu(change_in_SD_within_poor_environment))';
real_change_rich = nonzeros(triu(real_change_in_SD_within_rich))';
real_change_poor = nonzeros(triu(real_change_in_SD_within_poor))';
real_change_points(:,1,:) = SD_LT(:,2,:) - SD_LT(:,1,:); %mid-low
real_change_points(:,2,:) = SD_LT(:,3,:) - SD_LT(:,1,:); %high-low
real_change_points(:,3,:) = SD_LT(:,3,:) - SD_LT(:,2,:); %high-mid

plot(change_rich,'.--','Color',[0.7 0.3 0.3], 'LineWidth', 2, 'MarkerSize', 15); hold on
plot(change_poor,'.--','Color',[0.3 0.3 0.7], 'LineWidth', 2,'MarkerSize', 15);
plot(real_change_rich,'.-','Color',[0.7 0.3 0.3], 'LineWidth', 2,'MarkerSize', 15) % rich
plot(real_change_poor,'.-','Color',[0.3 0.3 0.7], 'LineWidth', 2,'MarkerSize', 15) % poor
scatter([1 2 3], real_change_points(:,:,1),[],[0.7 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.2) % rich
scatter([1 2 3], real_change_points(:,:,2),[],[0.3 0.3 0.7], 'filled', 'MarkerFaceAlpha', 0.2) % poor

set(gca,'XLim',[0 4],'XTick',[ 1 2 3],'XTickLabel',{'Mid-Low','High-low','High-mid'})
ylabel('Change in SD leaving time (s) (patch yield)')