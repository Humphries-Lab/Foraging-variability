%% ---- Wrapper script to simulate patch foraging task with stochastic models ---- %%
% Emma Scholey
% latest update 3 August 2023

clear
close all

addpath('./helperFunctions')
addpath('./analyticalComputation')

%% user options

% study options
study = 'leheron'; % study to simulate/fit data to. Options are leheron, contrerashuerta, kane

% model options
modelTable = readtable('./foragingModelTable.xlsx'); 
modelNum = 5; % model type - choose from foragingModelTable 

% simulation options
simOptions.type = 'simulate_new'; % 'simulate_new' if simulating from scratch, 'simulate_fit' if simulating pre-fit parameters for each subject
simOptions.blockPresentation = 'separate'; % either 'combined' (fit as one continuous task) or 'separate' (fit rich and poor environments as separate blocks)

% set parameters - options below will override if simulating already fit parameters
simOptions.nSim = 50;
simOptions.params.rich = [5, 0, 0.01, 1, 0, 0]; % {'beta', 'epsilon', 'alphaRho', 'alphaPatch', 'lambda', 'gamma'}
simOptions.params.poor = [5, 0, 0.01, 1, 0, 0]; % {'beta', 'epsilon', 'alphaRho', 'alphaPatch', 'lambda', 'gamma'}

%% set up

% load task
task = buildTask(study,simOptions.blockPresentation); % set up task structure 

% load dataframe container for simulations, according to 
allData = buildData(study,task,simOptions);

% load model
model = table2struct(modelTable(modelTable.modelNumber == modelNum,:));

% load agent parameters
allParams = buildParams(study,task,model,simOptions); clear params
model.paramNames = allParams.names;

simOptions.nSim = size(allData.nStates,1); % refresh nSim if loaded fit data

%% run simulations
modelLT_mean = zeros([task.nEnviron,task.nPatch,simOptions.nSim]); % store leaving times
modelLT_sd = zeros([task.nEnviron,task.nPatch,simOptions.nSim]); % store leaving times

for iS = 1:simOptions.nSim

    for iB = 1:task.nEnviron

        iS

        % what is the agent data for this block?
        agent.currentBlock = allData.blockOrder(iS,iB); % get current block for this subject based on their block order

        % set experienced avgRR depending on model
        switch model.rhoFunction
            case 'mvt'
                agent.experiencedAvgRR = task.optAvgRR; % MVT optimal
            case 'block'
                agent.experiencedAvgRR = allData.experiencedAvgRR(iS,:); % real experienced avgRR
        end
        
        agent.nStates = allData.nStates(iS,agent.currentBlock);
        agent.data = allData.data{agent.currentBlock}; % empty generic container
        agent.patchOrder = task.patchOrder{agent.currentBlock}'; % generic patch order for all participants
        agent.blockSwitchPatchN = 0; % only required for fitting 

        agentParams = table2array(allParams.params{agent.currentBlock}(iS,:)); % specific for the subject

        [NLL,results,out] = simulate_MVT_model(task,model,agent,agentParams,simOptions);

        % Extract leaving times (LT) for each patch and block type
        for iP = 1:task.nPatch
            modelLT_mean(agent.currentBlock,iP,iS) = mean(results.leaveT(results.patchOrder(1:end-1) == iP), 'omitnan');
            modelLT_sd(agent.currentBlock,iP,iS) = std(results.leaveT(results.patchOrder(1:end-1) == iP), 'omitnan');

        end

        % store timestep data 
        simData{iS,agent.currentBlock} = out;
    end
end

% save for paper figures
subject_leave_times_simulated.mean = modelLT_mean;
subject_leave_times_simulated.sd = modelLT_sd;

save_name = ['modelLT_', study,'_',simOptions.blockPresentation, '_M', sprintf('%d',modelNum), '.mat'];
save_path = '../data/simulation_data/';
save([save_path, save_name],'subject_leave_times_simulated');

%% plot against participant data

[dataLT] = loadLeaveT(study,task); % load mean leaving times for each subject x patch x env

% summarise model data
modelMean = mean(modelLT_mean,3);
subjectMean = mean(dataLT,3);

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
ylabel('Leaving time (s)')

legend({'model - rich', 'model - poor', 'data - rich', 'data - poor'})

% SD figure 
figure
plot(mean(modelLT_sd(1,:,:),3),'.-', 'Color',color.rich ,'LineWidth', 2, 'MarkerSize', 15); hold on
plot(mean(modelLT_sd(2,:,:),3),'.-', 'Color',color.poor,'LineWidth', 2, 'MarkerSize', 15)
set(gca,'XLim',[0 4],'XTick',1:task.nPatch,'XTickLabel',task.patchNames, 'YLim', [0 10])

%% plot avgRR and beta 

figure
plot(simData{6,1}.rho); ylim([0 30]); ylabel('avgRR'),xlabel('time in block'); hold on % rich
plot(simData{6,2}.rho); ylim([0 30]); ylabel('avgRR'),xlabel('time in block') % poor
legend({'rich', 'poor'})
figure
plot(simData{6,1}.beta); ylim([0 0.5]); ylabel('beta'),xlabel('time in block'); hold on % rich
plot(simData{6,2}.beta); ylim([0 0.5]); ylabel('beta'),xlabel('time in block') % poor
legend({'rich', 'poor'})

%% note to self 
% for CH and Kane datasets - more than one block of each environment 
% I fit as their actual block order e.g. rich/poor/rich/poor etc 
% but when I simulate, I join all rich together in a big block, and all
% poor together. This will need changing for simulating their results
% (don't expect it to make a big difference). See AET task. 
