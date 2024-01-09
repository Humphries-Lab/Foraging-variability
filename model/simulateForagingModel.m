%% ---- Wrapper script to simulate patch foraging task with stochastic models ---- %%
% Emma Scholey
% latest update 3 August 2023

clear
close all

addpath('./helperFunctions')
addpath('./analyticalComputation')

%% user options

% study options
simOptions.study = 'contrerashuerta'; % study to simulate/fit data to. Options are leheron, contrerashuerta, kane

% model options
modelTable = readtable('./foragingModelTable.xlsx'); 
modelNum = 27; % model type - choose from foragingModelTable 

% simulation options
simOptions.type = 'simulate_fit'; % 'simulate_new' if simulating from scratch, 'simulate_fit' if simulating pre-fit parameters for each subject
simOptions.blockPresentation = 'combined'; % either 'combined' (fit as one continuous task) or 'separate' (fit rich and poor environments as separate blocks)

% set parameters - options below will override if simulating already fit parameters
simOptions.nSim = 50;
simOptions.params.rich = [0.28, 0, 0.007, 1, 0, 0, -3]; % {'beta', 'epsilon', 'alphaRho', 'alphaPatch', 'lambda', 'gamma', 'bias'}
simOptions.params.poor = [0.31, 0, 0.0003, 1, 0, 0, -3]; % {'beta', 'epsilon', 'alphaRho', 'alphaPatch', 'lambda', 'gamma', 'bias'}

%% set up

% load task
task = buildTask(simOptions.study,simOptions.blockPresentation); % set up task structure 

% load dataframe container for simulations, according to 
allData = buildData(task,simOptions);

% load model
model = table2struct(modelTable(modelTable.modelNumber == modelNum,:));

% load agent parameters
allParams = buildParams(task,model,simOptions); clear params
model.paramNames = allParams.names;

%% run simulations
for iS = 1:allData.nSubj
    iS
    agent.experiencedAvgRR = allData.experiencedAvgRR(iS,:); % real experienced avgRR

    for iB = 1:task.nBlocks

        agent.nStates = allData.nStates;
        agent.data = allData.data; % empty generic container
        agent.currentEnv = allData.blockOrder(iS,iB);

        agent.patchOrder = task.patchOrder{agent.currentEnv}'; % generic patch order for all participants

        agentParams = table2array(allParams.params{allData.blockOrder(iS,iB)}(iS,:)); % specific for the subject, depending on environment
  
        [NLL,simLT{iS,iB},simData{iS,iB}] = simulate_MVT_model(task,model,agent,agentParams,simOptions);
    end
end

%% leaving times per environment x patch x subject
modelLT_mean = zeros([task.nEnviron,task.nPatch,allData.nSubj]); % store leaving times
modelLT_sd = zeros([task.nEnviron,task.nPatch,allData.nSubj]); % store leaving times

% Extract leaving times (LT) for each patch and block type
for iS = 1:allData.nSubj
    for iE = 1:task.nEnviron
        envData = vertcat(simLT{iS,allData.blockOrder(iS,:)==iE});
        for iP = 1:task.nPatch
            modelLT_mean(iE,iP,iS) = mean(envData.leaveT(envData.patchOrder(1:end-1) == iP), 'omitnan');
            modelLT_sd(iE,iP,iS) = std(envData.leaveT(envData.patchOrder(1:end-1) == iP), 'omitnan');
        end
    end
end

% save for paper figures
subject_leave_times_simulated.mean = modelLT_mean;
subject_leave_times_simulated.sd = modelLT_sd;

% save_name = ['modelLT_', simOptions.study,'_',simOptions.blockPresentation, '_M', sprintf('%d',modelNum), '.mat'];
% save_path = '../data/simulation_data/';
% save([save_path, save_name],'subject_leave_times_simulated');

%% plot against participant data

[dataLT, dataLT_SD] = loadLeaveT(simOptions.study,task); % load mean leaving times for each subject x patch x env

% summarise model data
modelMean = mean(modelLT_mean,3);
subjectMean = mean(dataLT,3);
subjectSD = mean(dataLT_SD,3);

figure

color.rich = [0.7 0.3 0.3];
color.poor = [0.3 0.3 0.7];

% plot the simulated data
plot(modelMean(1,:),'.--', 'Color',color.rich ,'LineWidth', 2, 'MarkerSize', 15), hold on
plot(modelMean(2,:),'.--', 'Color',color.poor,'LineWidth', 2, 'MarkerSize', 15)

% plot the actual experimental data
plot(subjectMean(1,:),'.-', 'Color',color.rich ,'LineWidth', 2, 'MarkerSize', 15)
plot(subjectMean(2,:),'.-', 'Color',color.poor,'LineWidth', 2, 'MarkerSize', 15)

set(gca,'XLim',[0 4],'XTick',1:task.nPatch,'XTickLabel',task.patchNames, 'YLim', [0 30])
ylabel('Mean leaving time (s)')
title('Model', modelNum)
legend({'model - rich', 'model - poor', 'data - rich', 'data - poor'})

% SD figure 
figure
plot(mean(modelLT_sd(1,:,:),3),'.--', 'Color',color.rich ,'LineWidth', 2, 'MarkerSize', 15); hold on
plot(mean(modelLT_sd(2,:,:),3),'.--', 'Color',color.poor,'LineWidth', 2, 'MarkerSize', 15)
plot(subjectSD(1,:),'.-', 'Color',color.rich ,'LineWidth', 2, 'MarkerSize', 15)
plot(subjectSD(2,:),'.-', 'Color',color.poor,'LineWidth', 2, 'MarkerSize', 15)
set(gca,'XLim',[0 4],'XTick',1:task.nPatch,'XTickLabel',task.patchNames, 'YLim', [0 10])
title('Model', modelNum)
ylabel('Mean SD of leaving time (s)')


%% plot avgRR and beta 
% REQUIRES UPDATING NOW THAT MORE BLOCKS HAVE BEEN ADDED FOR CH/Kane
% figure
% plot(simData{6,1}.rho); ylim([0 30]); ylabel('avgRR'),xlabel('time in block'); hold on % rich
% plot(simData{6,2}.rho); ylim([0 30]); ylabel('avgRR'),xlabel('time in block') % poor
% legend({'rich', 'poor'})
% figure
% plot(simData{6,1}.beta); ylim([0 0.5]); ylabel('beta'),xlabel('time in block'); hold on % rich
% plot(simData{6,2}.beta); ylim([0 0.5]); ylabel('beta'),xlabel('time in block') % poor
% legend({'rich', 'poor'})
