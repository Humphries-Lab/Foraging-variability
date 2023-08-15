function df = buildData(study,task,funcOptions)

% build dataframe, df
% funcOptions - simulate_new: simulate from scratch, simulate_fit: simulate best fit parameters; fit: use real data

% specify actions
leave = 1;
stay = 2;

switch funcOptions.type
    case 'simulate_new'

        % generate fake agent session
        df.blockOrder = [repmat([1 2], [funcOptions.nSim/2, 1]); repmat([2 1], [funcOptions.nSim/2, 1])];
        df.experiencedAvgRR = repmat([task.optAvgRR], [funcOptions.nSim,1]); % just assume experiencedAvgRR is the same as optimal, if simulating from scratch 
        df.nStates = repmat(task.blockTime, [funcOptions.nSim 2]);

        % set up empty container to store their data
        df.data = cell([1 2]);

        for iB = 1:task.nEnviron
            df.data{iB}.action = zeros(task.blockTime+1,1); % what action taken
            df.data{iB}.rho = zeros(task.blockTime+1,1); % estimated averageRR
            df.data{iB}.estPatchRR = zeros(task.blockTime+1,1); % estimated patchRR

            df.data{iB} = struct2table(df.data{iB});
        end

    case 'simulate_fit'

        % load the agent sessions for specific study
        switch study
            case 'leheron'
                load ../data/experiment_data/leheron_experiencedAvgRR.mat
                load ../data/experiment_data/leheron_blockorder.mat
                numSubjects = size(BlockOrder,1);
            case 'contrerashuerta'
                numSubjects = 29;
                % assume half of participants see rich vs poor block first
                BlockOrder = [repmat([1 2], [30/2, 1]); repmat([2 1], [30/2, 1])];

                experiencedAvgRR = readtable('../data/experiment_data/contrerashuerta_experiencedAvgRR.csv'); % load their real average RR in the 2 blocks
                df.experiencedAvgRR = [experiencedAvgRR.Rich_Self,experiencedAvgRR.Poor_Self];

            case 'kane'
                numSubjects = 8;
                % assume half of participants see rich vs poor block first
                BlockOrder = [repmat([1 2], [numSubjects/2, 1]); repmat([2 1], [numSubjects/2, 1])];

                T = readtable('../data/experiment_data/kane2019-rats-fig-1-data.csv');
                T = T(contains(T.Experiment,'Travel'),:); % only looking at travel time experiment

                numSubjects = size(unique(T.Subject),1);
                subID = unique(unique(T.Subject));

                for iS = 1:numSubjects % for each subject
                    subjLT = T(T.Subject == subID(iS),:); % extract their summarised leaving times
                    % load each subjects' experienced AvgRR for combined fitting
                    for iE = 1:task.nEnviron % for each environment, find their average reward rate
                        % take the average of averageRR across the 5
                        % testing days
                        currentEnv = subjLT(subjLT.Travel == task.travelTime(iE)*10,:);
                        endTrials = find(currentEnv.X == 1)-1; % find indices for the last trial of each session
                        endTrials = endTrials(2:end); % remove first trial
                        endTrials(end+1) = numel(currentEnv.X); % add on last trial

                        % find experienced avgRR for each day 
                        currentAvgRR = currentEnv.CumulativeReward_mL(endTrials)./currentEnv.CumulativeTime(endTrials);
                        % average over the 5 testing days 
                        experiencedAvgRR(iS,iE) = mean(currentAvgRR)*10000; 
                    end
                end

        end

        % set up empty container to store their data
        df.data = cell([1 2]);

        for iB = 1:task.nEnviron
            df.data{iB}.action = zeros(task.blockTime+1,1); % what action taken
            df.data{iB}.rho = zeros(task.blockTime+1,1); % estimated averageRR
            df.data{iB}.estPatchRR = zeros(task.blockTime+1,1); % estimated patchRR

            df.data{iB} = struct2table(df.data{iB});
        end

        df.nStates = repmat(task.blockTime, [numSubjects 2]);
        df.experiencedAvgRR = experiencedAvgRR;
        df.blockOrder = BlockOrder;

    case 'fit'

        switch study
            case 'leheron'
                % get subject data
                subj = [1:23 25:40];
                numSubjects = size(subj,2);

                load('../data/experiment_data/leheron_experiencedAvgRR.mat'); % load their real average RR in the 2 blocks
                df.experiencedAvgRR = experiencedAvgRR;

                load ../data/experiment_data/leheron_t_young.mat % load patch leaving times for each participant

                for iS = 1:numSubjects % for each subject

                    subjLT = t_young(t_young.subj == iS,:); % extract their summarised leaving times

                    for iB = 1:task.blockNum  % 1 if fitting together, 2 if fitting separately

                        if strcmp(funcOptions.blockPresentation,'separate')
                            leaveT = subjLT.leaveT(subjLT.env == iB); % pull out leaving times
                            patchOrder = subjLT.patch(subjLT.env == iB); % pull out patch order
                            env = subjLT.env(subjLT.env == iB); % pull out environment
                        elseif strcmp(funcOptions.blockPresentation,'combined')
                            leaveT = subjLT.leaveT; % pull out leaving times - note that this will do it in the correct block order for the participant [1 2] or [2 1]
                            patchOrder = subjLT.patch; % pull out patch order
                            env = subjLT.env; % pull out environment
                        end

                        leaveT = round(leaveT);

                        a = cell([numel(leaveT),1]);
                        % transform leaving times into stay/leave actions for each state
                        for ii = 1:numel(leaveT)
                            a{ii} = repelem([stay leave], [leaveT(ii) task.travelTime(env(ii))]);
                        end

                        A = cat(2, a{:})'; % concatenate all actions

                        df.nStates(iS,iB) = numel(A);
                        df.patchOrder{iS,iB} = patchOrder;
                        df.leaveT{iS,iB} = leaveT;
                        df.env{iS,iB} = env;
                        df.nObservations(iS,iB) = sum(A == 2); % only count stay states as an observation (they can't make choices whilst leaving)

                        df.data{iS}{iB}.action = [A;nan]; % what action taken
                        df.data{iS}{iB}.rho = zeros(numel(A)+1,1); % estimated averageRR
                        df.data{iS}{iB}.estPatchRR = zeros(numel(A)+1,1); % estimated patchRR

                        df.data{iS}{iB} = struct2table(df.data{iS}{iB});
                    end

                end

            case 'contrerashuerta'

                T = readtable('../data/experiment_data/contrerashuerta_trialbytrial_exp2.csv');
                T = T(T.ben == 1,:); % self condition only
                numSubjects = size(unique(T.sub),1);


                experiencedAvgRR = readtable('../data/experiment_data/contrerashuerta_experiencedAvgRR.csv'); % load their real average RR in the 2 blocks
                df.experiencedAvgRR = [experiencedAvgRR.Rich_Self,experiencedAvgRR.Poor_Self];

                for iS = 1:numSubjects % for each subject

                    subjLT = T(T.sub == iS,:); % extract their summarised leaving times

                    for iB = 1:task.blockNum  % 1 if fitting together, 2 if fitting separately

                        if strcmp(funcOptions.blockPresentation,'separate')
                            leaveT = subjLT.leaveT(subjLT.env == iB); % pull out leaving times
                            patchOrder = subjLT.patch(subjLT.env == iB); % pull out patch order
                            env = subjLT.env(subjLT.env == iB); % pull out environment
                        elseif strcmp(funcOptions.blockPresentation,'combined')
                            leaveT = subjLT.leaveT; % pull out leaving times - note that this will do it in the correct block order for the participant [1 2] or [2 1]
                            patchOrder = subjLT.patch; % pull out patch order
                            env = subjLT.env; % pull out environment
                        end

                        leaveT = round(leaveT);

                        a = cell([numel(leaveT),1]);
                        % transform leaving times into stay/leave actions for each state
                        for ii = 1:numel(leaveT)
                            a{ii} = repelem([stay leave], [leaveT(ii) task.travelTime(env(ii))]);
                        end

                        A = cat(2, a{:})'; % concatenate all actions

                        df.nStates(iS,iB) = numel(A);
                        df.patchOrder{iS,iB} = patchOrder;
                        df.leaveT{iS,iB} = leaveT;
                        df.env{iS,iB} = env;
                        df.nObservations(iS,iB) = sum(A == 2); % only count stay states as an observation (they can't make choices whilst leaving)

                        df.data{iS}{iB}.action = [A;nan]; % what action taken
                        df.data{iS}{iB}.rho = zeros(numel(A)+1,1); % estimated averageRR
                        df.data{iS}{iB}.estPatchRR = zeros(numel(A)+1,1); % estimated patchRR

                        df.data{iS}{iB} = struct2table(df.data{iS}{iB});
                    end

                end


            case 'kane'
                T = readtable('../data/experiment_data/kane2019-rats-fig-1-data.csv');
                T = T(contains(T.Experiment,'Travel'),:); % only looking at travel time experiment

                numSubjects = size(unique(T.Subject),1);
                subID = unique(unique(T.Subject));

                for iS = 1:numSubjects % for each subject

                    subjLT = T(T.Subject == subID(iS),:); % extract their summarised leaving times

                    for iB = 1:task.blockNum  % 1 if fitting together, 2 if fitting separately

                        % make dataset consistent to leheron and
                        % contreras-huerta
                        if strcmp(funcOptions.blockPresentation,'separate')

                            leaveT = subjLT.StateInPatch(subjLT.Travel == task.travelTime(iB)*10 & subjLT.Decision == 1); % pull out leaving times per patch
                            patchOrder = subjLT.startVolume(subjLT.Travel == task.travelTime(iB)*10 & subjLT.Decision == 1); % pull out patch order
                            patchOrder(patchOrder == task.r0(1)/1000) = 1; % convert patches to [1 2 3] 
                            patchOrder(patchOrder == task.r0(2)/1000) = 2;
                            patchOrder(patchOrder == task.r0(3)/1000) = 3;

                            env = repelem(iB,numel(leaveT))'; % dummy code environment [rich = 1, poor = 2]

                        elseif strcmp(funcOptions.blockPresentation,'combined')

                            leaveT = subjLT.StateInPatch(subjLT.Decision == 1); % pull out leaving times per patch
                            patchOrder = subjLT.startVolume(subjLT.Decision == 1); % pull out patch order
                            patchOrder(patchOrder == task.r0(1)/1000) = 1; % convert patches to [1 2 3] 
                            patchOrder(patchOrder == task.r0(2)/1000) = 2;
                            patchOrder(patchOrder == task.r0(3)/1000) = 3;
                            env = subjLT.Travel(subjLT.Decision == 1); % dummy code environment [rich = 1, poor = 2]
                            env(env == task.travelTime(1)*10) = 1;
                            env(env == task.travelTime(2)*10) = 2;

                        end

                        leaveT = round(leaveT);

                        a = cell([numel(leaveT),1]);
                        % transform leaving times into stay/leave actions for each state
                        for ii = 1:numel(leaveT)
                            a{ii} = repelem([stay leave], [leaveT(ii) task.travelTime(env(ii))]);
                        end

                        A = cat(2, a{:})'; % concatenate all actions

                        df.nStates(iS,iB) = numel(A);
                        df.patchOrder{iS,iB} = patchOrder;
                        df.leaveT{iS,iB} = leaveT;
                        df.env{iS,iB} = env;
                        df.nObservations(iS,iB) = sum(A == 2); % only count stay states as an observation (they can't make choices whilst leaving)

                        df.data{iS}{iB}.action = [A;nan]; % what action taken
                        df.data{iS}{iB}.rho = zeros(numel(A)+1,1); % estimated averageRR
                        df.data{iS}{iB}.estPatchRR = zeros(numel(A)+1,1); % estimated patchRR

                        df.data{iS}{iB} = struct2table(df.data{iS}{iB});
                   
                    end

                    for iE = 1:task.nEnviron % for each environment, find their average reward rate
                        % take the average of averageRR across the 5
                        % testing days
                        currentEnv = subjLT(subjLT.Travel == task.travelTime(iE)*10,:);
                        endTrials = find(currentEnv.X == 1)-1; % find indices for the last trial of each session
                        endTrials = endTrials(2:end); % remove first trial
                        endTrials(end+1) = numel(currentEnv.X); % add on last trial

                        % find experienced avgRR for each day 
                        currentAvgRR = currentEnv.CumulativeReward_mL(endTrials)./currentEnv.CumulativeTime(endTrials);
                        % average over the 5 testing days 
                        df.experiencedAvgRR(iS,iE) = mean(currentAvgRR)*10000; 
                    end                                   
                end


        end

end

