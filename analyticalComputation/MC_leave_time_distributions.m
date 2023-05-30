% script to Monte Carlo generate leave-time distributions for softmax
% policy
% Sweep over range of beta parameter, for 3 different patch types (r0)
% 
% Mark Humphries 13 May 2023. 

clearvars
close all;

number_samples = 250;  % number of Monte Carlo samples per parameter pair

beta = logspace(-3,0);  % space of softmax temperatures to calculate; 
                        % maximum of beta=1 here, as actual temperature is
                        % beta * r0
%epsilon = 0.01; 

t_max = 100;   % maximum time in patch (for explicit calculations)
t_step = 1;  % time-step at which to calculate estimates of E and VAR
alpha = 0.075;  % decay rate from LeHeron et al 202
r0 = [32.5 45 57.5];  % initial patch values from LeHeron et al 2020

n_steps = round(t_max / t_step);  % number of time-steps


%% MC sampling
t_series = (t_step* (1:n_steps))';  % the sequence of time-steps tested

sample_leave_time = zeros(numel(beta),numel(r0),number_samples);
E_leave = zeros(numel(beta),numel(r0)); VAR_leave = zeros(numel(beta),numel(r0));
for iB = 1:numel(beta)
    for iR = 1:numel(r0)
       %% calculate expected probability of leaving on each time-step
        f_leave = zeros(n_steps,1); p = zeros(n_steps,1); r = zeros(n_steps,1);
        for n = 1:n_steps
            time_now = n * t_step;                          % what is actual time on this time-step?          
            r(n) = reward_at_t_exp(time_now,r0(iR),alpha);  % reward on that time-step
            p(n) = p_leave_softmax(r(n),beta(iB));          % probability on that time-step
            %p(n) = p_leave_lapse(r(n),beta(iB), epsilon);          % probability on that time-step

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

        E_leave(iB,iR) = sum(t_series.*f_leave);
        VAR_leave(iB,iR) = sum((t_series - E_leave(iB,iR)).^2 .* f_leave); 
        
        %% Monte Carlo sampling of distribution
        for iSample = 1:number_samples
            blnLeave = 0; n = 1;
            while ~blnLeave
                time_now = n * t_step;                          % what is actual time on this time-step?          
                r = reward_at_t_exp(time_now,r0(iR),alpha);  % reward on that time-step
                p = p_leave_softmax(r,beta(iB));          % probability on that time-step
                 %p = p_leave_lapse(r,beta(iB), epsilon);          % probability on that time-step
                
                blnLeave = rand <= p;
                n = n + 1;
            end
            sample_leave_time(iB,iR,iSample) = time_now;
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
semilogx(beta,E_leave,'k'); hold on
semilogx(beta,E_MC_leave); hold on
xlabel('Beta: higher = exploit')
ylabel('Expected leaving time (s)')

figure
semilogx(beta,SD_leave,'k'); hold on
semilogx(beta,SD_MC_leave);
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
% 
% figure 
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

% E and SD at participant fit 
E_leave(ixPoor,:)
E_leave(ixRich,:)

SD_leave_mean(1,:) = SD_leave(ixRich,:);
SD_leave_mean(2,:) = SD_leave(ixPoor,:);

 %% compare to participant data
% 
% % initialise
% numBlocks = 2;
% numPatches = 3;
% 
% load t_young.mat % load summarised subject data 
% numSubjects = numel(unique(t_young.subj));
% 
% % test analytical predictions 
% 
% SD_LT = zeros([numSubjects, numPatches, numBlocks]); % each row is subject SD, each column patch type, 3d is block type
% 
% for iSubj = 1:numSubjects % for each subject
%     for iPatch = 1:numPatches % each patch
%         for iBlock = 1:numBlocks % each block
%             SD_LT(iSubj,iPatch,iBlock) = std(t_young.leaveT(t_young.subj == iSubj & t_young.env == iBlock & t_young.patch == iPatch));
%         end
%     end
% end
% 
% % between environments 
% 
% real_change_in_SD_between_environments = mean(SD_LT(:,:,2) - SD_LT(:,:,1)); % poor - rich
% % our results agree in that there is higher SD in poor environment compared
% % to rich, but differences in real SDs are much smaller than analytical SDs
% % could this be related to noisier SDs in participants, which mask any
% % differences?
% % and don't see increase in SD differences as patch yield increases
% % e.g. (0.5, 0.03 0.22 participants) compared to (4.5, 7, 9 analytical)  
% 
% % within environments 
% 
% real_change_in_SD_within_rich = mean(SD_LT(:,:,1)) - mean(SD_LT(:,:,1))';
% real_change_in_SD_within_poor = mean(SD_LT(:,:,2)) - mean(SD_LT(:,:,2))';
% 
% %mean SDs
% SD_LT_mean = squeeze(mean(SD_LT))';
% 
% % plot model vs participant SDs across patches and environment
% figure
% 
% plot(SD_leave(ixRich,:),'.--','Color',[0.7 0.3 0.3], 'LineWidth', 2, 'MarkerSize', 15); hold on
% plot(SD_leave(ixPoor,:),'.--','Color',[0.3 0.3 0.7], 'LineWidth', 2,'MarkerSize', 15);
% plot(SD_LT_mean(1,:),'.-','Color',[0.7 0.3 0.3], 'LineWidth', 2,'MarkerSize', 15) % rich
% plot(SD_LT_mean(2,:),'.-','Color',[0.3 0.3 0.7], 'LineWidth', 2,'MarkerSize', 15) % poor
% scatter([1 2 3], SD_LT(:,:,1),[],[0.7 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.1) % rich
% scatter([1 2 3], SD_LT(:,:,2),[],[0.3 0.3 0.7], 'filled', 'MarkerFaceAlpha', 0.1) % poor
% 
% set(gca,'XLim',[0 4],'XTick',[ 1 2 3],'XTickLabel',{'Low','Medium','High'})
% set(gca,'YLim', [0 12])
% ylabel('SD of leaving time (s)')
% 
% % plot between environment differences 
% figure
% 
% plot(change_in_SD_between_environments,'.--', 'Color', 'black', 'LineWidth', 2, 'MarkerSize', 15); hold on
% plot(real_change_in_SD_between_environments,'.-', 'Color', 'black', 'LineWidth', 2,'MarkerSize', 15) % rich
% scatter([1 2 3], SD_LT(:,:,2) - SD_LT(:,:,1), 'black',  'filled', 'MarkerFaceAlpha', 0.1) % rich
% 
% set(gca,'XLim',[0 4],'XTick',[ 1 2 3],'XTickLabel',{'Low','Medium','High'})
% ylabel('Change in SD leaving time (s) (poor - rich)')
% 
% 
% % plot within environment differences 
% figure
% change_rich = nonzeros(triu(change_in_SD_within_rich_environment))';
% change_poor = nonzeros(triu(change_in_SD_within_poor_environment))';
% real_change_rich = nonzeros(triu(real_change_in_SD_within_rich))';
% real_change_poor = nonzeros(triu(real_change_in_SD_within_poor))';
% real_change_points(:,1,:) = SD_LT(:,2,:) - SD_LT(:,1,:); %mid-low
% real_change_points(:,2,:) = SD_LT(:,3,:) - SD_LT(:,1,:); %high-low
% real_change_points(:,3,:) = SD_LT(:,3,:) - SD_LT(:,2,:); %high-mid
% 
% plot(change_rich,'.--','Color',[0.7 0.3 0.3], 'LineWidth', 2, 'MarkerSize', 15); hold on
% plot(change_poor,'.--','Color',[0.3 0.3 0.7], 'LineWidth', 2,'MarkerSize', 15);
% plot(real_change_rich,'.-','Color',[0.7 0.3 0.3], 'LineWidth', 2,'MarkerSize', 15) % rich
% plot(real_change_poor,'.-','Color',[0.3 0.3 0.7], 'LineWidth', 2,'MarkerSize', 15) % poor
% scatter([1 2 3], real_change_points(:,:,1),[],[0.7 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.2) % rich
% scatter([1 2 3], real_change_points(:,:,2),[],[0.3 0.3 0.7], 'filled', 'MarkerFaceAlpha', 0.2) % poor
% 
% set(gca,'XLim',[0 4],'XTick',[ 1 2 3],'XTickLabel',{'Mid-Low','High-low','High-mid'})
% ylabel('Change in SD leaving time (s) (patch yield)')