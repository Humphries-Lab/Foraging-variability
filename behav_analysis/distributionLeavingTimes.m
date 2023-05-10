% analysis of distribution of leaving times for Le Heron data 
% Emma Scholey 4 May 2023
clear
close all;
addpath(genpath('~/Dropbox/foraging/code'))

SetUpEnviron % get the environment parameters
blockNames = {'rich', 'poor'};

%% analysis of real subject data 
% get their data
subj = [1:23 25:40];

% initialise
numSubjects = size(subj,2);
numBlocks = size(blockNames, 2);

%subj = 1:4;
load ~/Dropbox/foraging/raw_data/summary/young_variables/t_young.mat
group = {'Young_HC_ReDo/yc%g_forage.mat', subj};

% calculate each subject's variances in leaving times in rich vs poor
% environments 

distrLT = zeros([numSubjects, 2]); % each row is subject variance, each column is rich vs poor block
for iSubj = 1:numSubjects % for each subject 
    subjLT_rich = t_young.leaveT(t_young.subj == iSubj & t_young.env == 1);
    subjLT_poor = t_young.leaveT(t_young.subj == iSubj & t_young.env == 2);
    distrLT(iSubj,1) = std(subjLT_rich);
    distrLT(iSubj,2) = std(subjLT_poor);
end

close all

% plot of distribution of leaving times - not normal (tail) 
figure
h1 = histogram(distrLT(:,1), 8);
hold on
h2 = histogram(distrLT(:,2), 8);
xlabel('Standard deviation of leaving times')
ylabel('Frequency of observations')
h1.FaceColor = [0.6353    0.0784    0.1843];
h2.FaceColor = [0    0.4471    0.7412];
legend({'Rich', 'Poor'})

distrLT_diff = distrLT(:,1) - distrLT(:,2); % differences between rich and poor block
% higher values indicate more variance in leaving times in rich environment
figure
histogram(distrLT_diff)
sum(distrLT_diff > 0) % 17 subjects had more variance in rich environment compared to poor environment 
sum(distrLT_diff < 0) % 22 subjects had more variance in poor environment compared to rich environment 
xlabel('Differences in std of leaving times')
ylabel('Frequency of observations')

signrank(distrLT(:,1), distrLT(:,2)); % no significant difference in std between rich and poor blocks
mean(distrLT(:,1))
mean(distrLT(:,2))

%% analysis of distributions of simulated data 

load ~/Dropbox/foraging/outputs/M2/SimulationResults.mat
distrLT = zeros([size(SimData, 1), 2]); % each row is subject variance, each column is rich vs poor block
for iSim = 1:size(SimData, 1) % for each simulation 
    simLT_rich = SimData{iSim}{1}.LeavingTime;
    simLT_poor = SimData{iSim}{2}.LeavingTime;
    distrLT(iSim,1) = std(simLT_rich);
    distrLT(iSim,2) = std(simLT_poor);
end

close all

% plot of distribution of leaving times - not normal (tail) 
figure
h1 = histogram(distrLT(:,1), 8);
hold on
h2 = histogram(distrLT(:,2), 8);
xlabel('Standard deviation of leaving times')
ylabel('Frequency of observations')
h1.FaceColor = [0.6353    0.0784    0.1843];
h2.FaceColor = [0    0.4471    0.7412];
legend({'Rich', 'Poor'})

distrLT_diff = distrLT(:,1) - distrLT(:,2); % differences between rich and poor block
% higher values indicate more variance in leaving times in rich environment
figure
histogram(distrLT_diff)
sum(distrLT_diff > 0) % 17 subjects had more variance in rich environment compared to poor environment 
sum(distrLT_diff < 0) % 22 subjects had more variance in poor environment compared to rich environment 
xlabel('Differences in std of leaving times')
ylabel('Frequency of observations')

[h p] = ttest(distrLT(:,1), distrLT(:,2)) % no significant difference in std between rich and poor blocks
mean(distrLT(:,1))
mean(distrLT(:,2))
