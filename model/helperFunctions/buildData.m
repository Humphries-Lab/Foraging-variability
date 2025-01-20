function df = buildData(task,funcOptions)

% Foraging variability project 
% Emma Scholey, 16 January 2025 

% BUILDDATA generates a container in the correct format for model fitting
% and simulation

% INPUTS:
%   task: provides the task variables specific to each dataset, such as 
%   optimal average reward rates and time in block
%
%   funcOptions: based on whether you are simulating (requiring an empty 
%   container), or fitting (prepares actual subject leaving times to be in 
%   the correct format for fitting per timestep rather than per patch.

% specify actions
leave = 1;
stay = 2;

% load real data to set up simulations
load(sprintf('../data/experiment_data/%s/%s_subj_var',funcOptions.study,funcOptions.study));

% get number of subjects to simulate, and their block orders/average reward rates
switch funcOptions.type
    case 'simulate_new'
        df.nSubj = funcOptions.nSim;
        df.blockOrder = subj_var.blockOrder(randi(size(subj_var.blockOrder,1),[df.nSubj 1]),:); % for each fake participant, create fake block order based on real participant data
        df.experiencedAvgRR = repmat([task.optAvgRR], [df.nSubj,1]); % just assume experiencedAvgRR is the same as optimal, if simulating from scratch

    case {'simulate_fit','fit'}
        df.nSubj = size(subj_var.blockOrder,1);
        df.blockOrder = subj_var.blockOrder; % take actual block order and experienced average RR of real subjects
        df.experiencedAvgRR = subj_var.experiencedAvgRR;
end

% build dataframe (for each second in the task)

switch funcOptions.type
    case {'simulate_new','simulate_fit'} % if simulating
        df.data.action = zeros(task.blockTime+1,1); % what action taken
        df.data.rho = zeros(task.blockTime+1,1); % estimated averageRR
        df.data.estPatchRR = zeros(task.blockTime+1,1); % estimated patchRR
        df.data = struct2table(df.data);

        df.nStates = task.blockTime; % number of timesteps in each block

    case 'fit'

        trialLeaveT = readtable(sprintf('../data/experiment_data/%s/%s_trialbytrial.csv',funcOptions.study,funcOptions.study));
        % find the patch number where the block switches

        if strcmp(funcOptions.study, 'contrerashuerta')
            trialLeaveT = trialLeaveT(trialLeaveT.ben == 1,:); % exclude other condition
        end

        for iS = 1:df.nSubj % for each subject

            subjTrialLeaveT = trialLeaveT(trialLeaveT.sub == iS,:); % extract their summarised leaving times

            leaveT = subjTrialLeaveT.leaveT; % pull out leaving times - note that this will do it in the correct block order for the participant
            patchOrder = subjTrialLeaveT.patch; % pull out patch order
            env = subjTrialLeaveT.env; % pull out environment
            switchIndex = subj_var.blockSwitchIndex{iS}; % index informing which patch to re-initialise estimates (i.e. when new block starts)

            leaveT = round(leaveT);
            a = cell([numel(leaveT),1]);

            % transform leaving times into stay/leave actions for each
            % timestep
            for ii = 1:numel(leaveT)
                a{ii} = repelem([stay leave], [leaveT(ii) task.travelTime(env(ii))]);
            end

            A = cat(2, a{:})'; % concatenate all actions

            df.nStates(iS) = numel(A);
            df.patchOrder{iS} = patchOrder;
            df.leaveT{iS} = leaveT;
            df.env{iS} = env;
            df.nObservations(iS) = sum(A == 2); % only count stay timesteps as an observation (they can't make choices whilst leaving)
            df.switchIndex{iS} = switchIndex;
            df.data{iS}.action = [A;nan]; % what action taken
            df.data{iS}.rho = zeros(numel(A)+1,1); % estimated averageRR
            df.data{iS}.estPatchRR = zeros(numel(A)+1,1); % estimated patchRR

            df.data{iS} = struct2table(df.data{iS});
        end

end

