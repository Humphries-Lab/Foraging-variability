%% ---- Wrapper script to simulate patch foraging task with stochastic models ---- %%
% Emma Scholey 9 Jun 2022
% latest update 16 January 2025


clear
close all

addpath('./helperFunctions')
addpath(genpath('../../data/experiment_data/'))

%% user options

nRun = 50; % how many iterations to simulate for each agent
modelNum = 17; % model number - choose from foragingModelTable 

% study options
simOptions.study = 'contrerashuerta'; % study to simulate/fit data to. Options are leheron (field-human), contrerashuerta (berry-human), or kane (rat)

% simulation options
simOptions.type = 'simulate_fit'; % 'simulate_new' if simulating from scratch, 'simulate_fit' if simulating pre-fit parameters for each subject

% set parameters - options below will override if simulating already fit parameters
simOptions.nSim = 39;
simOptions.params = [0 0 0 0.011705 0 31.021 -1.5274 -2.1147 0 0]; % {'beta','beta_rich', 'beta_poor', 'alphaRho', 'alphaPatch', 'lambda', 'gamma', 'bias', 'bias_rich','bias_poor'}

%% set up

% model options
modelTable = readtable('./foragingModelTable.xlsx'); 

% load task
task = buildTask(simOptions.study); % set up task structure 

% load dataframe container for simulations
allData = buildData(task,simOptions);

% load model
model = table2struct(modelTable(modelTable.modelNumber == modelNum,:));

% load agent parameters
allParams = buildParams(model,simOptions); clear params
model.paramNames = allParams.names;

%% run simulations

% empty containers
run_modelLT_mean = zeros(task.nEnviron, task.nPatch, allData.nSubj, nRun);
run_modelLT_SD = zeros(task.nEnviron, task.nPatch, allData.nSubj, nRun);

run_beta = zeros(task.blockTime, task.nEnviron, allData.nSubj, nRun);
run_rho = zeros(task.blockTime+1, task.nEnviron, allData.nSubj, nRun);

run_early_SD = zeros(allData.nSubj, 2, nRun);
run_late_SD = zeros(allData.nSubj, 2, nRun);

for iR = 1:nRun

    iR 
    for iS = 1:allData.nSubj

        agent.experiencedAvgRR = allData.experiencedAvgRR(iS,:); % real experienced avgRR (calculated post-hoc)

        initialiseRho = mean(task.optAvgRR); % set up rho for block 1 (only needed for learning models)

        for iB = 1:task.nBlocks % why are we doing by block? because we don't know how many trials (patches) in each time-restricted block, so we need to simulate each block separately

            agent.nStates = allData.nStates;
            agent.data = allData.data; % empty generic container
            agent.data.rho(1) = initialiseRho; % carry over rho from previous block

            agent.currentEnv = allData.blockOrder(iS,iB);

            agent.patchOrder = task.patchOrder{agent.currentEnv}'; % generic patch order for all participants

            agentParams = allParams.params{iS,:}; % specific for the subject, depending on environment

            [NLL,simLT{iS,iB},simData{iS,iB}] = simulate_MVT_model(task,model,agent,agentParams,simOptions); % run the simulation

            initialiseRho = simData{iS,iB}.rho(end); % what was the last rho from previous block (to carry over)
        end

    end
   

    %% leaving times per environment x patch x subject
    modelLT_mean = zeros([task.nEnviron,task.nPatch,allData.nSubj]); % store leaving times
    modelLT_SD = zeros([task.nEnviron,task.nPatch,allData.nSubj]); % store leaving times
    betas = zeros([task.blockTime,task.nEnviron,allData.nSubj]); % store leaving times
    rho = zeros([task.blockTime+1,task.nEnviron,allData.nSubj]); % store leaving times

    % Extract leaving times (LT) and trajectories for each patch and block type
    for iS = 1:allData.nSubj
        for iE = 1:task.nEnviron
            envData = vertcat(simLT{iS,allData.blockOrder(iS,:)==iE});
            for iP = 1:task.nPatch
                modelLT_mean(iE,iP,iS) = mean(envData.leaveT(envData.patchOrder(1:end-1) == iP), 'omitnan');
                modelLT_SD(iE,iP,iS) = std(envData.leaveT(envData.patchOrder(1:end-1) == iP), 'omitnan');
            end
        end

        data_rich = [simData{iS,allData.blockOrder(iS,:)==1}]; 
        data_poor = [simData{iS,allData.blockOrder(iS,:)==2}];

        betas_rich = mean([data_rich.beta],2);
        betas_poor = mean([data_poor.beta],2);

        rho_rich = mean([data_rich.rho],2);
        rho_poor = mean([data_poor.rho],2);

        betas(:,:,iS) = [betas_rich, betas_poor];
        rho(:,:,iS) = [rho_rich, rho_poor];

    end
        
    run_modelLT_mean(:,:,:,iR) = modelLT_mean;
    run_modelLT_SD(:,:,:,iR) = modelLT_SD;

    run_beta(:,:,:,iR) = betas;
    run_rho(:,:,:,iR) = rho;


end

% average across iterations 
simulated_leave_times.mean = mean(run_modelLT_mean,4);
simulated_leave_times.sd = mean(run_modelLT_SD,4);

trajectories.beta = mean(run_beta,4);
trajectories.rho = mean(run_rho,4);

%save the simulated data for each model
 % save_name = ['sim_LT_', simOptions.study, '_M', sprintf('%d',modelNum), '.mat'];
 % save_path = '../data/simulation_data/';
 % save([save_path, save_name],'simulated_leave_times', 'trajectories');

%% plot against participant data
    load(sprintf('../data/experiment_data/%s/subject_LT_%s', simOptions.study,simOptions.study), '-mat');

% summarise model data

group_model_mean = mean(simulated_leave_times.mean,3);
group_model_SD = mean(simulated_leave_times.sd,3);

% summarise subject data
group_subject_mean = mean(subject_leave_times.mean,3);
group_subject_SD = mean(subject_leave_times.sd,3);

figure

color.rich = [0.7 0.3 0.3];
color.poor = [0.3 0.3 0.7];

% plot the simulated data
plot(group_model_mean(1,:),'.--', 'Color',color.rich ,'LineWidth', 2, 'MarkerSize', 15), hold on
plot(group_model_mean(2,:),'.--', 'Color',color.poor,'LineWidth', 2, 'MarkerSize', 15)

% plot the actual experimental data
plot(group_subject_mean(1,:),'.-', 'Color',color.rich ,'LineWidth', 2, 'MarkerSize', 15)
plot(group_subject_mean(2,:),'.-', 'Color',color.poor,'LineWidth', 2, 'MarkerSize', 15)

set(gca,'XLim',[0 4],'XTick',1:task.nPatch,'XTickLabel',task.patchNames, 'YLim', [0 30])
ylabel('Mean leaving time (s)')
title('Model', modelNum)
legend({'model - rich', 'model - poor', 'data - rich', 'data - poor'})
set(findall(gcf,'-property','FontSize'),'FontSize',18)

% SD figure 
figure
plot(group_model_SD(1,:),'.--', 'Color',color.rich ,'LineWidth', 2, 'MarkerSize', 15); hold on
plot(group_model_SD(2,:),'.--', 'Color',color.poor,'LineWidth', 2, 'MarkerSize', 15)
plot(group_subject_SD(1,:),'.-', 'Color',color.rich ,'LineWidth', 2, 'MarkerSize', 15)
plot(group_subject_SD(2,:),'.-', 'Color',color.poor,'LineWidth', 2, 'MarkerSize', 15)
set(gca,'XLim',[0 4],'XTick',1:task.nPatch,'XTickLabel',task.patchNames, 'YLim', [0 10])
title('Model', modelNum)
ylabel('Mean SD of leaving time (s)')
set(findall(gcf,'-property','FontSize'),'FontSize',18)

