%% Script to plot data from foraging simulations (MVT or R-learning model)
% Emma Scholey

% clear up
clear all
load all_models.mat


%% run plots for all models
for model = 25
    % check if folder exists
    if not(isfolder(sprintf('~/Dropbox/foraging-project/results/M%d', model)))
        mkdir(sprintf('~/Dropbox/foraging-project/results/M%d', model))
    end

    %% plot individual plots of leaving stats
    save_plots = 0;
    load('~/Dropbox/foraging-project/data/raw_data/summary/young_variables/subMean_young.mat')

    plot_leaving_stats_line(all_models_gT_behaviour{model}, subMean_young)
    if save_plots == 1
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])
        print(sprintf('~/Dropbox/foraging-project/results/M%d/M%d_leaving_times', model, model), '-r600','-dpng')
    end

    %% calculate, store and plot Q stats
    save_plots = 0;
    if ~ismember(model, [1 2])
        % average and plot Q table across runs for each patch/environ
        plot_Q(all_models_gT_Q_stay{model}, all_models_gT_Q_leave{model}, model);
        if save_plots == 1
            print(sprintf('~/Dropbox/foraging-project/results/M%d/M%d_Q_values', model, model), '-dpng')
        end
    end
end

%% subplot leaving stats
% select which models to select from first, so not subplotting too many
subplot_leaving_stats_line(all_models_gT_behaviour);
if save_plots == 1
    print('~/Dropbox/foraging-project/results/subplot_leaving_rr', '-dpdf','-fillpage')
end

%% plot reward rate function for patches

RR_Leave = [21.8678 18.5632]; % this is the average background RR at which subjects should leave, for Rich and Poor env... -
... respectively,  derived from running Nils' script with 32.5 45 57.5 as A, 0.075 as a, and a 6 sec delay between patches
    ... this value also corresponds to the instantaneous RR on the gain function at which subjects should leave.
    clear optLT
A=[32.5 45 57.5];a=0.075;
for e=1:2 %each env
    for p=1:3 % each patch type
        optLT(e,p)=(log(RR_Leave(e)/A(p)))/-a;
    end
end

environ.time_step = 1; % adjust time step if required

for j=1:length(A) % for each patch type
    fun = @(T) A(j)*exp(-T*0.075); %Exponential function currently
    rr(:,j) = fun(0:50);
    k = 0:0.05:50;
    for i = 1:length(k)
        r(i,j) = integral(fun,k(i),k(i)+environ.time_step); % This creates matrix 'q' with each cell reward earned in time_step epoch at that reward rate
    end
end

c = [0 0 0;0.4  0.4  0.4; 0.7    0.7    0.7];
figure
plot(0:50, rr(:,1), 'LineWidth', 2, 'Color', c(1,:))
hold on
plot(0:50, rr(:,2), 'LineWidth', 2, 'Color', c(2,:))
plot(0:50, rr(:,3), 'LineWidth', 2, 'Color', c(3,:))
scatter(optLT(2,:), RR_Leave(2), [], [0    0.4471    0.7412], 'filled')
scatter(optLT(1,:), RR_Leave(1), [], [0.6353    0.0784    0.1843], 'filled')
hold off

xlabel('Time in patch (s)', 'FontSize', 14)
ylabel('Foreground reward rate', 'FontSize', 14)
poor_ref = refline(0, environ.averageRR.poor);
rich_ref = refline(0, environ.averageRR.rich);

poor_ref.Color = [0    0.4471    0.7412];
rich_ref.Color = [0.6353    0.0784    0.1843];

poor_ref.LineStyle = '--';
rich_ref.LineStyle = '--';
title('Reward rate of different patch types', 'FontSize', 16)
legend('Low yield', 'Medium yield', 'High yield')

% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])
% print('/Users/exs165/Dropbox/foraging-project/results/rr_function', '-dsvg')


%%
% % inferential stats on experiment data
% poor = [subMean_poor, ones([length(subMean_poor), 1])];
% rich = [subMean_rich, 2*ones([length(subMean_rich), 1])];
%     y = array2table([poor; rich]);
%     tmp = stack(y, 1:3);
%     g1 = tmp.Var4;
%     g2 = tmp.Var1_Var2_Var3_Indicator;
%     [p, tbl, stats] = anovan(tmp.Var1_Var2_Var3,{g1, g2}, 'model', 'interaction');

%% plot which model gets the highest reward

% extract total reward earned in Le Heron - code below from his script
group = {'/Users/exs165/Dropbox/foraging-project/data/raw_data/Young_HC_ReDo/yc%g_forage.mat', [1:23 25:40]};
alldata = {};
nn=length(group{2})
for k = 1:length(group{2}) %work through all subjects in group
    data = load(sprintf(group{1},group{2}(k))); %load subject data
    d = data.result.data; %retrive trial data
    %can remove data fields here if wish to
    alldata{k} = d;
end

for k = 1:length(group{2})
    total_reward_rich(k,:) = sum([alldata{k}([alldata{k}.block] == 1).patchReward]);
    total_reward_poor(k,:) = sum([alldata{k}([alldata{k}.block] == 2).patchReward]);
end

ppt_reward_rich = mean(total_reward_rich);
ppt_reward_poor = mean(total_reward_poor);

mean_sum_reward = cat(1, all_models_mean_reward{:});

% add in confidence intervals

% plot poor environment
close all
figure
fig1 = bar(3:26, mean_sum_reward(3:26,2));
fig1.FaceColor = [0    0.4471    0.7412];
ylim([0, 21000])
hline = yline(mean_sum_reward(1,2), '--','Optimal reward - poor environment');
hline2 = yline(ppt_reward_poor, '-','Participant reward - poor environment');
hline.Color = [0    0.4471    0.7412];
hline.LineWidth = 2;
hline2.LineWidth = 2;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])
xlabel('Model number')
ylabel('Total reward obtained (ml)')
title('Poor environment')
print('~/Dropbox/foraging-project/results/reward_obtained_poor', '-dpng')

% plot rich environment
figure
fig2 = bar(3:26, mean_sum_reward(3:26,1));
fig2.FaceColor = [0.6353    0.0784    0.1843];
ylim([0, 21000])
hline = yline(mean_sum_reward(1,1), '--','Optimal reward - rich environment');
hline2 = yline(ppt_reward_rich, '-','Participant reward - rich environment');
hline.Color = [0.6353    0.0784    0.1843];
hline.LineWidth = 2;
hline2.LineWidth = 2;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])
xlabel('Model number')
ylabel('Total reward obtained (ml)')
title('Rich environment')
print('~/Dropbox/foraging-project/results/reward_obtained_rich', '-dpng')

%% To do

% - which model learns the average reward rate.
% - replicate main task exactly e.g. include practice set of 10 patches in
% each environment, and have rich then poor or poor then rich, to see if any
% block learning effects e.g. assume average rr continues over blocks