%% figure - plot total reward
%% THIS NEEDS FINISHING 

% plot which model gets the highest reward
clear 
model = 25;

%% load subject data 

% extract total reward earned in Le Heron - code below from his script
group = {'~/Dropbox/foraging/raw_data/Young_HC_ReDo/yc%g_forage.mat', [1:23 25:40]};
alldata = {};
nn=length(group{2})
for k = 1:length(group{2}) %work through all subjects in group
    data = load(sprintf(group{1},group{2}(k))); %load subject data
    d = data.result.data; %retrive trial data
    %can remove data fields here if wish to
    alldata{k} = d;
end

for k = 1:length(group{2})
    for i = 1:2 % [rich, poor]
    SubTotalReward(k,i) = sum([alldata{k}([alldata{k}.block] == i).patchReward]);
    end
end

SubMeanReward = mean(SubTotalReward);

%% load model reward
load(sprintf('~/Dropbox/foraging/outputs/M%d/SimulationResults.mat', model)) % load simulated results for each model

% plot poor environment
close all
figure
fig1 = bar(model, RewardMean(3:26,2));
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
hline2 = yline(SubMeanReward, '-','Participant reward - rich environment');
hline.Color = [0.6353    0.0784    0.1843];
hline.LineWidth = 2;
hline2.LineWidth = 2;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])
xlabel('Model number')
ylabel('Total reward obtained (ml)')
title('Rich environment')
print('~/Dropbox/foraging-project/results/reward_obtained_rich', '-dpng')