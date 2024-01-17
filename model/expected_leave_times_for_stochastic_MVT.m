% script to generate leave-time distributions for stochastic choice models
%
% Sweep over range of exploration parameter, for different patch types (task.r0)
% 
% 6 August 2023 Emma Scholey: include contreras-huerta task and load real
% task data for comparison
% 31 July 2023 Mark Humphries: efficient code (removed loop over n time-steps)
% 28 July 2023 Mark Humphries: original

clearvars
close all;

addpath('helperFunctions')

% specify study
study = 'kane'; 

% parameters of choice model
model = 'softmax'; % 'softmax', 'e-greedy', 'lapse'
explore_parameter = logspace(-3,0,100);  % space of softmax temperatures to calculate; 
                        % maximum of beta=1 here, as actual temperature is
                        % beta * task.r0         
%lapse_parameter = 0.04; % required for lapse model only 
bias_parameter = -2.5; % required for softmax model - intercept for over/under-harvesting
t_max = 1000;   % maximum time in patch (for explicit calculations)

switch study
    case 'leheron'
        task  = buildTask(study,'separate');
    case 'contrerashuerta'
        task  = buildTask(study,'separate');
    case 'kane'
        task  = buildTask(study,'separate');
end

color.rich = [0.10,0.62,0.47];
color.poor = [0.46,0.44,0.70];

% parameters of E and VAR calculation
t_step = 1;  % time-step at which to calculate estimates of E and VAR; 1 = trials (discrete Bernoulli process); <1 = approximation to continuous time

%% computes E and VAR as function of exploration parameter
n_steps = round(t_max / t_step);  % number of time-steps
t_series = (t_step* (1:n_steps))';  % the sequence of time-steps tested

E_leave = zeros(numel(explore_parameter),numel(task.r0)); VAR_leave = zeros(numel(explore_parameter),numel(task.r0)); f_leave_all = zeros(numel(explore_parameter),numel(task.r0),n_steps);

for iR = 1:numel(task.r0)
    % calculate reward function for this task.r0
    switch task.rewardFunction
        case 'exponential'
            reward_ts = reward_at_t_exp(t_series,task.r0(iR),task.decayRate);
        case 'linear'
            reward_ts = reward_at_t_linear(t_series,task.r0(iR),task.decayRate); 
        otherwise
            error('Unrecognised reward function')
    end

    % calculate expected probability of leaving on each time-step
    for iB = 1:numel(explore_parameter)
        switch model
            case 'softmax'
                p_leave_at_n = p_leave_softmax(reward_ts,explore_parameter(iB), bias_parameter,0);
            case 'e-greedy'
                p_leave_at_n = explore_parameter(iB) * ones(size(reward_ts)); 
            case 'lapse'
                p_leave_at_n = p_leave_lapse(reward_ts,explore_parameter(iB), lapse_parameter);
            otherwise
                error('Unrecognised model')
        end
        
        % correct probability of leaving by time-step
        p_leave_at_n = p_leave_at_n .* t_step;

        % probability of staying up to time-step n
        cumulative_p_of_staying_until_n = cumprod(1-p_leave_at_n);  

        % probability of staying up to time-step n-1 and leaving on time-step n
        f_leave = [p_leave_at_n(1); cumulative_p_of_staying_until_n(1:end-1).*p_leave_at_n(2:end)];  
        
        % store
        f_leave_all(iB,iR,:) = f_leave;
        
        % calculate E and VAR
        E_leave(iB,iR) = sum(t_series.*f_leave);
        VAR_leave(iB,iR) = sum((t_series - E_leave(iB,iR)).^2 .* f_leave); 
        
    end
end

%% derived statistics from calculations
SD_leave = sqrt(VAR_leave);
CV_leave = E_leave ./ SD_leave;

%% save data for paper figures 
data.E_leave = E_leave;
data.SD_leave = SD_leave;
data.beta = explore_parameter;

% save_name = ['expectedLT_', model, '_',study,'_04', '.mat'];
% save_path = '../../data/analytical_data/';
% save([save_path, save_name],'data');

%% plot results
E_fig = figure;
semilogx(explore_parameter,E_leave); hold on
switch model
    case {'softmax','lapse'}
        xlabel('Beta (higher = exploit)')
    case 'e-greedy'
        xlabel('\epsilon (higher = explore)')
        %set(gca,'XDir','reverse')
end

ylabel('Expected leaving time (s)')

SD_fig = figure;
semilogx(explore_parameter,SD_leave); hold on
switch model
    case {'softmax','lapse'}
        xlabel('Beta (higher = exploit)')
    case 'e-greedy'
        xlabel('\epsilon (higher = explore)')
        %set(gca,'XDir','reverse')
end
ylabel('SD of leaving time (s)')

%% extract expected beta values for specific study

[subjectLT_mean, subjectLT_sd, subjectRR_mean] = loadLeaveT(study,task); % load mean and std of leaving times for each subject x patch x env

% save for paper figures
subject_leave_times.mean = subjectLT_mean;
subject_leave_times.sd = subjectLT_sd;

% save_name = ['subjectLT_', study, '.mat'];
% save_path = '../../data/analytical_data/';
% save([save_path, save_name],'subject_leave_times');

% group averages 
mean_leave_data = mean(subjectLT_mean,3); 
SD_leave_data = mean(subjectLT_sd,3); 

% find corresponding values of beta for expected leaving time of mid patch 
ixRich = find(E_leave(:,2) >= mean_leave_data(1,2),1,"first");
ixPoor = find(E_leave(:,2) >= mean_leave_data(2,2),1,"first");

beta_rich = explore_parameter(ixRich);
beta_poor = explore_parameter(ixPoor);

%% plot predictions of means and SDs for these beta values 

% draw as coloured lines on above E and SD figs
figure(E_fig)
line([beta_rich beta_rich],[0 E_fig.Children.YLim(2)],'Color',color.rich)
line([beta_poor beta_poor],[0 E_fig.Children.YLim(2)],'Color',color.poor)

% plot lines on SD too...
figure(SD_fig)
line([beta_rich beta_rich],[0 E_fig.Children.YLim(2)],'Color',color.rich)
line([beta_poor beta_poor],[0 E_fig.Children.YLim(2)],'Color',color.poor)
ylim([0 10])

%% plot expected leaving time for these beta values, compared to experimental data
figure
plot(1:numel(task.r0),E_leave(ixRich,:),'.--','Color',color.rich,'LineWidth',1.5); hold on
plot(1:numel(task.r0),E_leave(ixPoor,:),'.--','Color',color.poor,'LineWidth',1.5)

% plot against experimental data
% plot(1:numel(task.r0),mean_leave_data(1,:),'.-','Color',color.rich,'LineWidth',1.5); hold on
% plot(1:numel(task.r0),mean_leave_data(2,:),'.-','Color',color.poor,'LineWidth',1.5)

% plot against optimal MVT 
plot(1:numel(task.r0),mean_leave_data(1,:),'.-','Color',color.rich,'LineWidth',1.5); hold on
plot(1:numel(task.r0),mean_leave_data(2,:),'.-','Color',color.poor,'LineWidth',1.5)

ylabel('Expected leaving time (s)')
title(sprintf('beta rich = %f, beta poor = %f', round(beta_rich,3),round(beta_poor,3)))

set(gca,'XTick',1:numel(task.r0),'XTickLabel',task.patchNames,'XLim',[0 4],'YLim',[0 inf])


%% compute SD predictions for model, compared to participant data

% expected SD 
model_change_in_SD_between_environments = SD_leave(ixPoor,:) - SD_leave(ixRich,:); 
model_change_in_SD_within_rich = nonzeros(triu(SD_leave(ixRich,:) - SD_leave(ixRich,:)'))';
model_change_in_SD_within_poor = nonzeros(triu(SD_leave(ixPoor,:) - SD_leave(ixPoor,:)'))';

% compare to experiment data
real_change_in_SD_between_environments = SD_leave_data(2,:) - SD_leave_data(1,:); % poor - rich
real_change_in_SD_within_rich = nonzeros(triu(SD_leave_data(1,:) - SD_leave_data(1,:)'))';
real_change_in_SD_within_poor = nonzeros(triu(SD_leave_data(2,:) - SD_leave_data(2,:)'))';

%% plot model vs participant SDs
figure

% model 
plot(SD_leave(ixRich,:),'.--','Color',color.rich, 'LineWidth', 1.5, 'MarkerSize', 10); hold on
plot(SD_leave(ixPoor,:),'.--','Color',color.poor, 'LineWidth', 1.5,'MarkerSize', 10);

% experiment
plot(SD_leave_data(1,:),'.-','Color',color.rich, 'LineWidth', 1.5,'MarkerSize', 10)
plot(SD_leave_data(2,:),'.-','Color',color.poor, 'LineWidth', 1.5,'MarkerSize', 10)

% include distribution of experiment SDs
scatter(1:task.nPatch, squeeze(subjectLT_sd(1,:,:)),[],color.rich, 'filled', 'MarkerFaceAlpha', 0.2)
scatter(1:task.nPatch, squeeze(subjectLT_sd(2,:,:)),[],color.poor, 'filled', 'MarkerFaceAlpha', 0.2)

set(gca,'XLim',[0 4],'XTick',1:task.nPatch,'XTickLabel',task.patchNames)
ylabel('SD of leaving time (s)')

% %% plot between environment SD effects 
% figure
% 
% plot(model_change_in_SD_between_environments,'.--', 'Color', 'black', 'LineWidth', 2, 'MarkerSize', 10); hold on
% plot(real_change_in_SD_between_environments,'.-', 'Color', 'black', 'LineWidth', 2,'MarkerSize', 10)
% scatter(1:task.nPatch, squeeze(subjectLT_sd(2,:,:)-subjectLT_sd(1,:,:)),'k', 'filled', 'MarkerFaceAlpha', 0.2)
% 
% set(gca,'XLim',[0 4],'XTick',1:task.nPatch,'XTickLabel',task.patchNames)
% ylabel('Change in SD leaving time (s) (poor - rich)')
% %legend({'Model' 'Data'})
% 
% %% plot within environment SD effects 
% figure
% 
% for iE = 1:task.nEnviron
%     switch study
%         case 'contrerashuerta'
%             points_real_change_in_SD_within(:,1,iE) = squeeze(subjectLT_sd(iE,2,:)-subjectLT_sd(iE,1,:)); % high-low
%             patchDiffNames = {'High-low'};
%             nPatchDiff = 1;
%         otherwise
%             points_real_change_in_SD_within(:,1,iE) = squeeze(subjectLT_sd(iE,2,:)-subjectLT_sd(iE,1,:)); % mid-low
%             points_real_change_in_SD_within(:,2,iE) = squeeze(subjectLT_sd(iE,3,:)-subjectLT_sd(iE,1,:)); % high-low
%             points_real_change_in_SD_within(:,3,iE) = squeeze(subjectLT_sd(iE,3,:)-subjectLT_sd(iE,2,:)); % high-mid
%             patchDiffNames = {'Mid-Low','High-low','High-mid'};
%             nPatchDiff = 3;
%     end
% end
% 
% plot(model_change_in_SD_within_rich,'.--','Color',color.rich, 'LineWidth', 1.5, 'MarkerSize', 10); hold on
% plot(model_change_in_SD_within_poor,'.--','Color',color.poor, 'LineWidth', 1.5,'MarkerSize', 10);
% plot(real_change_in_SD_within_rich,'.-','Color',color.rich, 'LineWidth', 1.5,'MarkerSize', 10) % rich
% plot(real_change_in_SD_within_poor,'.-','Color',color.poor, 'LineWidth', 1.5,'MarkerSize', 10) % poor
% scatter(1:nPatchDiff, points_real_change_in_SD_within(:,:,1),[],color.rich, 'filled', 'MarkerFaceAlpha', 0.2)
% scatter(1:nPatchDiff, points_real_change_in_SD_within(:,:,2),[],color.poor, 'filled', 'MarkerFaceAlpha', 0.2)
% 
% set(gca,'XLim',[0 4],'XTick',1:nPatchDiff,'XTickLabel',patchDiffNames)
% ylabel('Change in SD leaving time (s) (patch yield)')
% %legend({'Model - rich', 'Model - poor', 'Data - rich', 'Data - poor'})
% 
% %% reward rate predictions - plot reward rates for range of explore parameter
% for iR = 1:numel(task.r0)
%     switch task.rewardFunction
%         case 'exponential'
%             RR_leave(:,iR) = reward_at_t_exp(E_leave(:,iR),task.r0(iR),task.decayRate);
%         case 'linear'
%             RR_leave(:,iR) = reward_at_t_linear(E_leave(:,iR),task.r0(iR),task.decayRate); 
%         otherwise
%             error('Unrecognised reward function')
%     end
% end
% 
% RR_fig = figure;
% semilogx(explore_parameter,RR_leave); hold on
% switch model
%     case {'softmax','lapse'}
%         xlabel('Beta (higher = exploit)')
%     case 'e-greedy'
%         xlabel('\epsilon (higher = explore)')
%         %set(gca,'XDir','reverse')
% end
% ylabel('Reward rate at expected leaving time (s)')
% 
% line([beta_rich beta_rich],[0 RR_fig.Children.YLim(2)],'Color',color.rich)
% line([beta_poor beta_poor],[0 RR_fig.Children.YLim(2)],'Color',color.poor)
% legend({'low','medium','high'})
% 
% % plot Kane's dataset reward rates 
% 
% figure 
% 
% mean_leave_RR_data = mean(subjectRR_mean,3); 
% plot(1:numel(task.r0),mean_leave_RR_data(1,:),'.-','Color',color.rich,'LineWidth',1.5); hold on
% plot(1:numel(task.r0),mean_leave_RR_data(2,:),'.-','Color',color.poor,'LineWidth',1.5)
% 
% ylabel('Reward rate at expected leaving time (ml)')
% 
% set(gca,'XTick',1:numel(task.r0),'XTickLabel',task.patchNames,'XLim',[0 4],'YLim',[0 0.12])
