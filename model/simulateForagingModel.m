%% ---- Script to simulate patch foraging task ---- %%
% Emma Scholey
% latest update 3 August 2023

clear
close all
addpath('.\Main\Sub')

%% user options

% study options
study = 'kane'; % study to simulate/fit data to. TO DO: If contreras-huerta or kane, then only checked working with model 1 so far

% model options
modelNum = 1; % model type - see model table to check number to choose
modelTable = readtable('/Users/exs165/Dropbox/foraging/code/foragingModelTable.xlsx'); % change directory

% simulation options
simOptions.type = 'simulate_new'; % 'simulate_new' if simulating new parameters, 'simulate_fit' if simulating already fit parameters for each subject
simOptions.blockPresentation = 'separate'; % either 'combined' (fit as one continuous task) or 'separate' (fit rich and poor as separate blocks)

% options below will override if simulating already fit parameters
simOptions.nSim = 50;
simOptions.params.rich = [0.04, 0.2, 0.0006, 1, 3, 0.2]; % {'beta', 'epsilon', 'alphaRho', 'alphaPatch', 'lambda', 'gamma'}
simOptions.params.poor = [0.06, 0.1, 0.0009, 1, 0.2, 0.2]; % {'beta', 'epsilon', 'alphaRho', 'alphaPatch', 'lambda', 'gamma'}

%% set up

% load task
task = buildTask(study,simOptions.blockPresentation); % set up task

% load dataframe container for simulations
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

        agent.experiencedAvgRR = allData.experiencedAvgRR(iS,:); % specific for the subject
        agent.nStates = allData.nStates(iS,agent.currentBlock);
        agent.data = allData.data{agent.currentBlock}; % empty generic container
        agent.patchOrder = task.patchOrder{agent.currentBlock}'; % generic patch order for all participants

        agentParams = table2array(allParams.params{agent.currentBlock}(iS,:)); % specific for the subject

        if iB == 2 && strcmp(simOptions.blockPresentation, 'combined') % specify start values if it's the 2nd block and we're doing combined fitting of both blocks together
            agent.data.rho(1) = out.rho(end);
        end

        [NLL,results,out] = simulate_MVT_model(task,model,agent,agentParams,simOptions);

        % Extract leaving times (LT) for each patch and block type
        for iP = 1:task.nPatch
            modelLT_mean(agent.currentBlock,iP,iS) = mean(results.leaveT(results.patchOrder == iP), 'omitnan');
            modelLT_sd(agent.currentBlock,iP,iS) = std(results.leaveT(results.patchOrder == iP), 'omitnan');

        end

    end
end

% save for paper figures
subject_leave_times_simulated.mean = modelLT_mean;
subject_leave_times_simulated.sd = modelLT_sd;

save_name = ['modelLT_', study, '_M', sprintf('%d',modelNum), '.mat'];
save_path = '/Users/exs165/Dropbox/foraging/paper/figs/data/';
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

