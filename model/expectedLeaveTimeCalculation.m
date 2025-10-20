%% script to generate leave-time distributions for stochastic choice models
%
% Sweep over range of exploration parameter, for different patch types (task.r0)

% 11 June 2025 Emma Scholey: change bounds for bias (c) predictions
% 16 January 2025 Emma Scholey: to remove unused lapse model
% 13 Feb 2024 Emma Scholey: add bias parameter predictions
% 6 August 2023 Emma Scholey: include contreras-huerta task
% 31 July 2023 Mark Humphries: efficient code (removed loop over n time-steps)
% 28 July 2023 Mark Humphries: original

clearvars
close all;

addpath('../model/helperFunctions')

% specify study
study = 'contrerashuerta'; % either 'leheron', 'contrerashuerta', or 'kane'

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
bias_parameter = -2.05; % can also set to mean bias term from model fits. These are [-2.22 Le Heron, -2.52 CH, -2.05 Kane]
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

save_name = ['expectedLT_',study,'_rangebeta_SD.mat'];
save_path = '../data/analytical_data/';
save([save_path, save_name],'data');

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

beta_parameter = 0.4;
%bias_parameter = -[logspace(2,0,100)]; %when beta = 0.4
%bias_parameter = [logspace(-1,1,100)];%when beta = 0
bias_parameter = [-60:0.01:3]
%bias_parameter = [0:0.01:3];
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

save_name = ['expectedLT_',study,'_rangebias.mat'];
save_path = '../data/analytical_data/';
save([save_path, save_name],'data');

%% plot results

signedlog = @(x) sign(x) .* log10(1 + abs(x));
bias_signedlog = signedlog(bias_parameter);

E_fig = figure;
%semilogx(bias_signedlog,E_leave); hold on
plot(bias_signedlog, E_leave); hold on
xlabel('Bias (higher = exploit)')
ylabel('Expected leaving time (s)')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
set(gca,'box','off')
xticks = [-60 -10 -1 0 1 3];  % Choose meaningful bias values
set(gca, 'XTick', signedlog(xticks))
set(gca, 'XTickLabel', arrayfun(@num2str, xticks, 'UniformOutput', false))
%
% SD_fig = figure;
% semilogx(bias_parameter,SD_leave); hold on
%
% xlabel('Bias (higher = exploit)')
% ylabel('SD of leaving time (s)')
% set(findall(gcf,'-property','FontSize'),'FontSize',18)
% set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
% set(gca,'box','off')




%% GENERATE DATA FOR SD HEATMAP FIGURES

% parameters of VAR calculation
t_step = 1;  
beta_parameter = logspace(-2,0,500);
bias_parameter = -6:0.1:1;
t_max = 1000;   % maximum time in patch (for explicit calculations)

    %% computes VAR as function of both exploration parameters
    n_steps = round(t_max / t_step);  % number of time-steps
    t_series = (t_step* (1:n_steps))';  % the sequence of time-steps tested

    E_leave = zeros(numel(bias_parameter),numel(beta_parameter),numel(task.r0));
    VAR_leave = zeros(numel(bias_parameter),numel(beta_parameter),numel(task.r0));
    f_leave_all = zeros(numel(bias_parameter),numel(beta_parameter), numel(task.r0),n_steps);

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
        for iBias = 1:numel(bias_parameter)

            for iBeta = 1:numel(beta_parameter)

                p_leave_at_n = p_leave_softmax(reward_ts,beta_parameter(iBeta), bias_parameter(iBias),0);

                % correct probability of leaving by time-step
                p_leave_at_n = p_leave_at_n .* t_step;

                % probability of staying up to time-step n
                cumulative_p_of_staying_until_n = cumprod(1-p_leave_at_n);

                % probability of staying up to time-step n-1 and leaving on time-step n
                f_leave = [p_leave_at_n(1); cumulative_p_of_staying_until_n(1:end-1).*p_leave_at_n(2:end)];

                % store
                f_leave_all(iBias,iBeta,iR,:) = f_leave;

                % calculate VAR
                E_leave(iBias,iBeta,iR) = sum(t_series.*f_leave);
                VAR_leave(iBias,iBeta,iR) = sum((t_series - E_leave(iBias,iBeta,iR)).^2 .* f_leave);

            end
        end
    end

    SD_leave = sqrt(VAR_leave);

    if strcmp(study, 'contrerashuerta')
        SD_diff = SD_leave(:,:,2) - SD_leave(:,:,1);
    else
        SD_diff = SD_leave(:,:,3) - SD_leave(:,:,1);
    end


%% save data for paper figures
data.SD_leave_heatmap = SD_leave;
data.SD_leave_heatmap_diff = SD_diff;
data.beta = beta_parameter;
data.bias = bias_parameter;
save_name = ['expected_SD_',study,'.mat'];
save_path = '../data/analytical_data/';
save([save_path, save_name],'data');
