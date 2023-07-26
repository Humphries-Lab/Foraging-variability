% script to Monte Carlo generate leave-time distributions for different
% action selection policies - softmax and epsilon-softmax
% Sweep over range of lambda parameter, for 3 different patch types (r0)
% 
% Mark Humphries 13 May 2023. 
% Emma Scholey updated - latest 11 July 2023

clearvars
close all;

number_samples = 250;  % number of Monte Carlo samples per parameter pair

lambda = logspace(-1,1.2,100);  % space of lambda values to calculate - based on range of participant fits
                        
t_max = 100;   % maximum time in patch (for explicit calculations)
t_step = 0.1;  % time-step at which to calculate estimates of E and VAR
alpha = 0.075;  % decay rate from LeHeron et al 202
r0 = [32.5 45 57.5];  % initial patch values from LeHeron et al 2020

n_steps = round(t_max / t_step);  % number of time-steps

%rho = [21.8678, 18.5632]; % optimal average RR (rho) as defined by MVT

load ~/Dropbox/foraging/raw_data/experiencedAvgRR.mat % mean experienced participant average RR
rho = mean(experiencedAvgRR); clear experiencedAvgRR

%% MC sampling
t_series = (t_step* (1:n_steps))';  % the sequence of time-steps tested

sample_leave_time = zeros(numel(lambda),numel(r0),number_samples);

for iL = 1:numel(lambda)
    for iR = 1:numel(r0)
       %% calculate expected probability of leaving on each time-step
        f_leave_rich = zeros(n_steps,1); p_rich = zeros(n_steps,1); r = zeros(n_steps,1);
        f_leave_poor = zeros(n_steps,1); p_poor = zeros(n_steps,1);
        for n = 1:n_steps
            time_now = n * t_step;                          % what is actual time on this time-step?          
            r(n) = reward_at_t_exp(time_now,r0(iR),alpha);  % reward on that time-step
            beta_rich = lambda(iL)*1/(rho(1)); % beta is some function of averageRR, based on rich environment
            beta_poor = lambda(iL)*1/(rho(2)); % beta is some function of averageRR, based on poor environment

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

        %% Monte Carlo sampling of distribution
%         for iSample = 1:number_samples
%             blnLeave = 0; n = 1;
%             while ~blnLeave
%                 time_now = n * t_step;                          % what is actual time on this time-step?          
%                 r = reward_at_t_exp(time_now,r0(iR),alpha);  % reward on that time-step
%                 beta = lambda(iL)*1/optRho(1); % beta is some function of averageRR, based on rich environment currently
%                 p = p_leave_softmax(r,beta);          % probability on that time-step
% 
%                 blnLeave = rand <= p;
%                 n = n + 1;
%             end
%             sample_leave_time(iL,iR,iSample) = time_now;
%         end
    end
end

%% derived statistics from calculations
SD_leave_rich = sqrt(VAR_leave_rich);
SD_leave_poor = sqrt(VAR_leave_poor);

% CV_leave = E_leave ./ SD_leave;

%% statistics from MC simulations
% E_MC_leave = mean(sample_leave_time,3);
% SD_MC_leave = std(sample_leave_time,[],3); 

%% plot results
figure
semilogx(lambda,E_leave_poor,'k'); hold on
%semilogx(lambda,E_MC_leave); hold on
xlabel('Lambda: higher = exploit')
ylabel('Expected leaving time (s)')

figure
semilogx(lambda,SD_leave_poor,'k'); hold on
%semilogx(lambda,SD_MC_leave);
xlabel('Lambda: higher = exploit')
ylabel('SD of leaving time (s)')

% figure
% semilogx(lambda,CV_leave)
% xlabel('Lambda: higher = exploit')
% ylabel('Coeff Var of leaving time (s)')

% plot examples like LeHeron task performance
% roughly from Fig 2A:
% rich: [11 16 21]
% poor: [15 20 25]

% find index that agrees with middle-valued patch (as it occured at same
% rate in both environments)
%ixRich = find(E_leave_rich(:,2) > 9.5,1,'first');  % for optimal LTs, not participant LTs
ixRich = find(E_leave_rich(:,2) > 16,1,'first');  % for participant LTs

lambda_rich = lambda(ixRich)

%ixPoor = find(E_leave_poor(:,2) > 11,1,'first'); % for optimal LTs, not participant LTs
ixPoor = find(E_leave_poor(:,2) > 19,1,'first'); % for participant LTs

lambda_poor = lambda(ixPoor);

lambda_rich == lambda_poor % these should be the same value - if not then somethings gone wrong

figure 
plot(E_leave_rich(ixRich,:),'.-','Color',[0.7 0.3 0.3]); hold on
plot(E_leave_poor(ixPoor,:),'.-','Color',[0.3 0.3 0.7]);

% overlay MVT optimal times 
RR_Leave = [21.8678 18.5632]; % this is the average background RR at which subjects should leave, for Rich and Poor env... - 
... respectively.
    
clear optLT
A=[32.5 45 57.5];a=0.075;
for e=1:2 %each env
    for p=1:3 % each patch type
        optLT(e,p)=(log(RR_Leave(e)/A(p)))/-a;
    end
end
plot(optLT(1,:),'--','Color',[0.7 0.3 0.3])
plot(optLT(2,:),'--','Color',[0.3 0.3 0.7])

% plot rough participant values for comparison

load('~/Dropbox/foraging/raw_data/summary/young_variables/subMean_young.mat') % load real data

% plot the actual experimental data
NSub = size(subMean_young, 3);
SubLT(:,:,1) = squeeze(subMean_young(1,:,:)); % rich environment
SubLT(:,:,2) = squeeze(subMean_young(2,:,:)); % poor environment
SubLTMean = [mean(SubLT(:,:,1), 2), mean(SubLT(:,:,2), 2)];
SubSEM = [std(SubLT(:,:,1),[],2)./sqrt(NSub), std(SubLT(:,:,2), [], 2)./sqrt(NSub)]; % Standard Error

ts = tinv([0.025  0.975],NSub-1);      % T-Score
SubCI = ts(2).*SubSEM; % CIs

plot(SubLTMean(:,1),':', 'Color',[0.7 0.3 0.3]);
hold on
plot(SubLTMean(:,2),':', 'Color',[0.3 0.3 0.7]);

legend({'Model - rich', 'Model - poor', 'MVT - rich', 'MVT - poor', 'Subj - rich', 'Subj - poor'})
set(gca,'XLim',[0 4],'XTick',[ 1 2 3],'XTickLabel',{'Low','Medium','High'}, 'YLim', [0 30])
ylabel('Expected leaving time (s)')

%% question: which is bigger - change in SD of patch-type between environments, 
% or between patch types in same environment?

% do for LeHeron estimates here - this is a more general question 

change_in_SD_between_environments = SD_leave_poor(ixPoor,:) - SD_leave_poor(ixRich,:)
change_in_SD_within_rich_environment = SD_leave_rich(ixRich,:) - SD_leave_rich(ixRich,:)'
change_in_SD_within_poor_environment = SD_leave_poor(ixPoor,:) - SD_leave_poor(ixPoor,:)'


SD_leave_mean(1,:) = SD_leave_rich(ixRich,:);
SD_leave_mean(2,:) = SD_leave_poor(ixPoor,:);

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
 % our results agree in that there is higher SD in poor environment compared
 % to rich, but differences in real SDs are much smaller than analytical SDs
 % could this be related to noisier SDs in participants, which mask any
 % differences?
 % and don't see increase in SD differences as patch yield increases
 % e.g. (0.5, 0.03 0.22 participants) compared to (4.5, 7, 9 analytical)

 % within environments

 real_change_in_SD_within_rich = mean(SD_LT(:,:,1)) - mean(SD_LT(:,:,1))';
 real_change_in_SD_within_poor = mean(SD_LT(:,:,2)) - mean(SD_LT(:,:,2))';

 %mean SDs
 SD_LT_mean = squeeze(mean(SD_LT))';

 % plot model vs participant SDs across patches and environment
 figure

 plot(SD_leave_rich(ixRich,:),'.--','Color',[0.7 0.3 0.3], 'LineWidth', 2, 'MarkerSize', 15); hold on
 plot(SD_leave_poor(ixPoor,:),'.--','Color',[0.3 0.3 0.7], 'LineWidth', 2,'MarkerSize', 15);
 plot(SD_LT_mean(1,:),'.-','Color',[0.7 0.3 0.3], 'LineWidth', 2,'MarkerSize', 15) % rich
 plot(SD_LT_mean(2,:),'.-','Color',[0.3 0.3 0.7], 'LineWidth', 2,'MarkerSize', 15) % poor
 scatter([1 2 3], SD_LT(:,:,1),[],[0.7 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.1) % rich
 scatter([1 2 3], SD_LT(:,:,2),[],[0.3 0.3 0.7], 'filled', 'MarkerFaceAlpha', 0.1) % poor

 set(gca,'XLim',[0 4],'XTick',[ 1 2 3],'XTickLabel',{'Low','Medium','High'})
 set(gca,'YLim', [0 12])
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
