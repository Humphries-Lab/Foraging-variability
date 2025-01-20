%% script to generate leave-time distributions for stochastic choice models
%
% Sweep over range of exploration parameter, for different patch types (task.r0)

% 16 January 2025 Emma Scholey: to remove unused lapse model
% 13 Feb 2024 Emma Scholey: add bias parameter predictions
% 6 August 2023 Emma Scholey: include contreras-huerta task
% 31 July 2023 Mark Humphries: efficient code (removed loop over n time-steps)
% 28 July 2023 Mark Humphries: original

clearvars
close all;

addpath('helperFunctions')

% specify study
study = 'leheron'; % either 'leheron', 'contrerashuerta', or 'kane'

switch study
    case 'leheron'
        task  = buildTask(study);
    case 'contrerashuerta'
        task  = buildTask(study);
    case 'kane'
        task  = buildTask(study);
end

color.rich = [0.10,0.62,0.47];
color.poor = [0.46,0.44,0.70];

% parameters of E and VAR calculation
t_step = 1;  % time-step at which to calculate estimates of E and VAR; 1 = trials (discrete Bernoulli process); <1 = approximation to continuous time

% parameters of choice model
model = 'softmax'; % 'softmax', 'e-greedy'

%% ---------------- Expected leave times for BETA parameter -----------------------------%%

beta_parameter = logspace(-3,1,1000);
bias_parameter = -2; % can also set to mean bias term from model fits. These are [-2.28 Le Heron, -2.55 CH, -2.04 Kane]
t_max = 1000;   % maximum time in patch (for explicit calculations)

%% computes E and VAR as function of exploration parameter
n_steps = round(t_max / t_step);  % number of time-steps
t_series = (t_step* (1:n_steps))';  % the sequence of time-steps tested

E_leave = zeros(numel(beta_parameter),numel(task.r0)); VAR_leave = zeros(numel(beta_parameter),numel(task.r0)); f_leave_all = zeros(numel(beta_parameter),numel(task.r0),n_steps);

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
    for iB = 1:numel(beta_parameter)
        switch model
            case 'softmax'
                p_leave_at_n = p_leave_softmax(reward_ts,beta_parameter(iB), bias_parameter,0);
            case 'e-greedy'
                p_leave_at_n = beta_parameter(iB) * ones(size(reward_ts));
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
data.explore = beta_parameter;

% save_name = ['expectedLT_',study,'_rangebeta_','bias=',num2str(bias_parameter), '.mat'];
% save_path = '../data/analytical_data/';
% save([save_path, save_name],'data');

%% plot results
E_fig = figure;
semilogx(beta_parameter,E_leave); hold on
% for quick inspect of betas
ixRich = find(data.E_leave >= 4.77,1,"first");
ixPoor = find(data.E_leave >= 6.56,1,"first");
beta_rich = data.explore(ixRich);
beta_poor = data.explore(ixPoor);
line([beta_rich beta_rich],[0 E_fig.Children.YLim(2)], 'Color', color.rich, 'LineStyle', '--')
line([beta_poor beta_poor],[0 E_fig.Children.YLim(2)], 'Color', color.poor, 'LineStyle', '--')
legend({'low','mid','high','beta rich', 'beta poor'})

switch model
    case {'softmax'}
        xlabel('Beta (higher = exploit)')
    case 'e-greedy'
        xlabel('\epsilon (higher = explore)')
        %set(gca,'XDir','reverse')
end

ylabel('Expected leaving time (s)')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
set(gca,'box','off')
xlim([beta_parameter(1),beta_parameter(end)])

SD_fig = figure;
semilogx(beta_parameter,SD_leave); hold on
line([beta_rich beta_rich],[0 SD_fig.Children.YLim(2)], 'Color', color.rich, 'LineStyle', '--')
line([beta_poor beta_poor],[0 SD_fig.Children.YLim(2)], 'Color', color.poor, 'LineStyle', '--')
legend({'low','mid','high','beta rich', 'beta poor'})
xlim([beta_parameter(1),beta_parameter(end)])

switch model
    case 'softmax'
        xlabel('Beta (higher = exploit)')
    case 'e-greedy'
        xlabel('\epsilon (higher = explore)')
        %set(gca,'XDir','reverse')
end
ylabel('SD of leaving time (s)')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
set(gca,'box','off')


%% ---------------- Expected leave times for BIAS parameter -----------------------------%%

beta_parameter = 0;
bias_parameter = [logspace(-1,1,100)]; %when beta = 0

t_max = 1000;   % maximum time in patch (for explicit calculations)

%% computes E and VAR as function of exploration parameter
n_steps = round(t_max / t_step);  % number of time-steps
t_series = (t_step* (1:n_steps))';  % the sequence of time-steps tested

E_leave = zeros(numel(bias_parameter),numel(task.r0)); VAR_leave = zeros(numel(bias_parameter),numel(task.r0)); f_leave_all = zeros(numel(bias_parameter),numel(task.r0),n_steps);

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
    for iB = 1:numel(bias_parameter)
        p_leave_at_n = p_leave_softmax(reward_ts,beta_parameter, bias_parameter(iB),0);

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
data.explore = bias_parameter;

% save_name = ['expectedLT_',study,'_rangebias_','beta=',num2str(beta_parameter), '.mat'];
% save_path = '../data/analytical_data/';
% save([save_path, save_name],'data');

%% plot results
E_fig = figure;
semilogx(bias_parameter,E_leave); hold on

xlabel('Bias (higher = exploit)')
ylabel('Expected leaving time (s)')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
set(gca,'box','off')

SD_fig = figure;
semilogx(bias_parameter,SD_leave); hold on

xlabel('Bias (higher = exploit)')
ylabel('SD of leaving time (s)')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
set(gca,'box','off')
