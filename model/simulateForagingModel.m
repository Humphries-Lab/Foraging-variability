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
modelNum = 7; % model type - choose from foragingModelTable 

% simulation options
simOptions.type = 'simulate_fit'; % 'simulate_new' if simulating from scratch, 'simulate_fit' if simulating pre-fit parameters for each subject

% set parameters - options below will override if simulating already fit parameters
simOptions.nSim = 50;
simOptions.params = [0.4 0.41 0.45 0 0 0 4 0 -1 0 0]; % {'beta','beta_rich', 'beta_poor','epsilon', 'alphaRho', 'alphaPatch', 'lambda', 'gamma', 'bias', 'bias_rich','bias_poor'}

%% set up

% load task
task = buildTask(simOptions.study); % set up task structure 

% load dataframe container for simulations, according to 
allData = buildData(task,simOptions);

% load model
model = table2struct(modelTable(modelTable.modelNumber == modelNum,:));

% load agent parameters
allParams = buildParams(task,model,simOptions); clear params
model.paramNames = allParams.names;

%% run simulations
for iS = 1%:allData.nSubj
    iS
    agent.experiencedAvgRR = allData.experiencedAvgRR(iS,:); % real experienced avgRR

    for iB = 1:task.nBlocks

        agent.nStates = allData.nStates;
        agent.data = allData.data; % empty generic container
        agent.currentEnv = allData.blockOrder(iS,iB);

        agent.patchOrder = task.patchOrder{agent.currentEnv}'; % generic patch order for all participants

        agentParams = allParams.params{iS,:}; % specific for the subject, depending on environment
  
        [NLL,simLT{iS,iB},simData{iS,iB}] = simulate_MVT_model(task,model,agent,agentParams,simOptions);
    end
end

%% leaving times per environment x patch x subject
modelLT_mean = zeros([task.nEnviron,task.nPatch,allData.nSubj]); % store leaving times
modelLT_sd = zeros([task.nEnviron,task.nPatch,allData.nSubj]); % store leaving times
modelLT_mean_all = zeros([1,allData.nSubj]); % store leaving times
modelLT_sd_all = zeros([1,allData.nSubj]); % store leaving times
modelLT_mean_env = zeros([task.nEnviron,allData.nSubj]); % store leaving times
modelLT_sd_env = zeros([task.nEnviron,allData.nSubj]); % store leaving times

% Extract leaving times (LT) for each patch and block type
for iS = 1%:allData.nSubj
    subjData = vertcat(simLT{iS,:});
    modelLT_mean_all(iS) = mean(subjData.leaveT);
    modelLT_sd_all(iS) = std(subjData.leaveT);

    betas(:,:,iS) = [simData{iS,1}.beta,simData{iS,2}.beta];

    for iE = 1:task.nEnviron
        envData = vertcat(simLT{iS,allData.blockOrder(iS,:)==iE});
        modelLT_mean_env(iE,iS) = mean(envData.leaveT);
        modelLT_sd_env(iE,iS) = std(envData.leaveT);
        for iP = 1:task.nPatch
            modelLT_mean(iE,iP,iS) = mean(envData.leaveT(envData.patchOrder(1:end-1) == iP), 'omitnan');
            modelLT_sd(iE,iP,iS) = std(envData.leaveT(envData.patchOrder(1:end-1) == iP), 'omitnan');
        end
    end
end


% %save for paper figures
% simulated_leave_times.mean = modelLT_mean;
% simulated_leave_times.sd = modelLT_sd;
% 
% save_name = ['modelLT_', simOptions.study, '_M', sprintf('%d',modelNum), '.mat'];
% save_path = '../data/simulation_data/';
% save([save_path, save_name],'simulated_leave_times');
% 
% % save for paper figures
% simulated_leave_times.mean = modelLT_mean_all;
% simulated_leave_times.sd = modelLT_sd_all;
% 
% save_name = ['modelLT_all_', simOptions.study, '_M', sprintf('%d',modelNum), '.mat'];
% save_path = '../data/simulation_data/';
% save([save_path, save_name],'simulated_leave_times');
% 
% % save for paper figures
% simulated_leave_times.mean = modelLT_mean_env;
% simulated_leave_times.sd = modelLT_sd_env;
% 
% save_name = ['modelLT_env_', simOptions.study, '_M', sprintf('%d',modelNum), '.mat'];
% save_path = '../data/simulation_data/';
% save([save_path, save_name],'simulated_leave_times');

% save the computed beta's for paper figures 
save_name = ['model_betas_', simOptions.study, '_M', sprintf('%d',modelNum), '.mat'];
save_path = '../data/simulation_data/';
save([save_path, save_name],'betas');


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
set(findall(gcf,'-property','FontSize'),'FontSize',18)

% SD figure 
figure
plot(mean(modelLT_sd(1,:,:),3),'.--', 'Color',color.rich ,'LineWidth', 2, 'MarkerSize', 15); hold on
plot(mean(modelLT_sd(2,:,:),3),'.--', 'Color',color.poor,'LineWidth', 2, 'MarkerSize', 15)
plot(subjectSD(1,:),'.-', 'Color',color.rich ,'LineWidth', 2, 'MarkerSize', 15)
plot(subjectSD(2,:),'.-', 'Color',color.poor,'LineWidth', 2, 'MarkerSize', 15)
set(gca,'XLim',[0 4],'XTick',1:task.nPatch,'XTickLabel',task.patchNames, 'YLim', [0 10])
title('Model', modelNum)
ylabel('Mean SD of leaving time (s)')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
