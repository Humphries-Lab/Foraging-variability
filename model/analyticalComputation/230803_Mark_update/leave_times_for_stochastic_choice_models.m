% script to generate leave-time distributions for stochastic choice models
%
% Sweep over range of exploration parameter, for different patch types (r0)
% 
% 31 July 2023: efficient code (removed loop over n time-steps)
% Mark Humphries 28 July 2023

clearvars
close all;

% exporting plots
export_path = '/Users/exs165/Dropbox/foraging/code/analyticalComputation/230803_Mark_update';

figsize = [20 20 7 7];
color.rich = [0.7 0.7 0.7];
color.poor = [0.8 0.7 0.5];

% parameters of choice model
% model = 'softmax'; % 'e-greedy', 
% explore_parameter = logspace(-3,0);  % space of softmax temperatures to calculate; 
%                         % maximum of beta=1 here, as actual temperature is
%                         % beta * r0            
% t_max = 100;   % maximum time in patch (for explicit calculations)
                        
 model = 'e-greedy';
 explore_parameter = logspace(-2,0);   % fixed probability of leaving 
 t_max = 1000;   % maximum time in patch (for explicit calculations)

% task parameters: LeHeron et al 2020
 reward_function = 'exponential';
 alpha = 0.075;  % decay rate from LeHeron et al 2020
 r0 = [32.5 45 57.5];  % initial patch values from LeHeron et al 2020


% task parameters: rat foraging study of Kane et al 2019 eLife
% reward_function = 'linear';
% alpha = 8;  % muL
% r0 = [60, 90, 120];  % muL

% parameters of E and VAR calculation
t_step = 1;  % time-step at which to calculate estimates of E and VAR; 1 = trials (discrete Bernoulli process); <1 = approximation to continuous time


%% computes E and VAR as function of exploration parameter
n_steps = round(t_max / t_step);  % number of time-steps
t_series = (t_step* (1:n_steps))';  % the sequence of time-steps tested

E_leave = zeros(numel(explore_parameter),numel(r0)); VAR_leave = zeros(numel(explore_parameter),numel(r0)); f_leave_all = zeros(numel(explore_parameter),numel(r0),n_steps);

for iR = 1:numel(r0)
    % calculate reward function for this r0
    switch reward_function
        case 'exponential'
            reward_ts = reward_at_t_exp(t_series,r0(iR),alpha);
        case 'linear'
            reward_ts = reward_at_t_linear(t_series,r0(iR),alpha); 
        otherwise
            error('Unrecognised reward function')
    end

    % calculate expected probability of leaving on each time-step
    for iB = 1:numel(explore_parameter)
        switch model
            case 'softmax'
                p_leave_at_n = p_leave_softmax(reward_ts,explore_parameter(iB));
            case 'e-greedy'
                p_leave_at_n = explore_parameter(iB) * ones(size(reward_ts)); 
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

%% plot results
E_fig = figure;
semilogx(explore_parameter,E_leave); hold on
switch model
    case 'softmax'
        xlabel('Beta (higher = exploit)')
    case 'e-greedy'
        xlabel('\epsilon (higher = explore)')
        %set(gca,'XDir','reverse')
end
ylabel('Expected leaving time (s)')
exportPPTfig(gcf,['E_leave_' model '_' reward_function],export_path)

SD_fig = figure;
semilogx(explore_parameter,SD_leave); hold on
switch model
    case 'softmax'
        xlabel('Beta (higher = exploit)')
    case 'e-greedy'
        xlabel('\epsilon (higher = explore)')
        %set(gca,'XDir','reverse')
end
ylabel('SD of leaving time (s)')
exportPPTfig(gcf,['SD_leave_' model '_' reward_function],export_path)

%% extract examples of rich/poor environments as lower and higher beta

beta_rich = 0.1;
beta_poor = 0.2;

% draw those as coloured lines on above E and SD figs
figure(E_fig)
line([beta_rich beta_rich],[0 E_fig.Children.YLim(2)],'Color',color.rich)
line([beta_poor beta_poor],[0 E_fig.Children.YLim(2)],'Color',color.poor)
exportPPTfig(gcf,['E_leave_lines_' model '_' reward_function],export_path)

% then plot those values of E across patches = show patch and environment
% effects!
ixRich = find(explore_parameter >= beta_rich,1,"first");
ixPoor = find(explore_parameter >= beta_poor,1,"first");

figure
plot(1:3,E_leave(ixRich,:),'.-','Color',color.rich); hold on
plot(1:3,E_leave(ixPoor,:),'.-','Color',color.poor)
ylabel('Expected leaving time (s)')
set(gca,'XTick',1:3,'XTickLabel',{'Low','Mid','High'},'XLim',[0 4],'YLim',[0 inf])

exportPPTfig(gcf,['Patch_and_environment_effect_' model '_' reward_function],export_path)

% plot lines on SD too...
figure(SD_fig)
line([beta_rich beta_rich],[0 E_fig.Children.YLim(2)],'Color',color.rich)
line([beta_poor beta_poor],[0 E_fig.Children.YLim(2)],'Color',color.poor)

exportPPTfig(gcf,['SD_leave_lines_' model '_' reward_function],export_path)

%% plot all elements of the model

if strcmp(model,'softmax')
    % reward function
    figure
    plot(t_series,reward_ts)
    xlabel('Trial (n)')
    ylabel('Reward')
    exportPPTfig(gcf,['example_reward_fcn_' reward_function],export_path)
    
    % softmax
    r_example = 0:50;
    figure
    plot(r_example,p_leave_softmax(r_example,0.2),'k'); hold on
    plot(r_example,p_leave_softmax(r_example,1),'Color',[0.7 0.7 0.7]);
    xlabel('Reward')
    ylabel('P(leave)')
    exportPPTfig(gcf,'example_softmax',export_path)

    % p(leave|n)
    figure
    plot(t_series,p_leave_softmax(reward_ts,0.2),'k'); hold on
    plot(t_series,p_leave_softmax(reward_ts,1),'Color',[0.7 0.7 0.7]);
    xlabel('Trial (n)')
    ylabel('P(leave|n)')
    exportPPTfig(gcf,['example_p_leave_at_n_' reward_function],export_path)
end

%% plot all p(leave = n) as beta increases, for one patch type

example_patch = 1;

figure
plot(t_series,squeeze(f_leave_all(:,example_patch,:)));
xlabel('Time')
ylabel('p(leave = n)')
title(['p(leave = n) for patch ' num2str(example_patch)])


