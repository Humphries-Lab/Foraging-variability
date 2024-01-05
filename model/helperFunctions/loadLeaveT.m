%% load participant data and summarise into leaving times
function [subLT_mean, subLT_sd, subRR_mean] = loadLeaveT(study,task)

switch study

    case {'leheron', 'contrerashuerta'}
        trialLeaveT = readtable(sprintf('../data/experiment_data/%s_trialbytrial.csv',study));

        if strcmp(study, 'contrerashuerta')
            trialLeaveT = trialLeaveT(trialLeaveT.ben == 1,:); % exclude other condition
        end

        nSub = numel(unique(trialLeaveT.sub));

        subLT_mean = zeros([task.nEnviron,task.nPatch,nSub]);
        subLT_sd = zeros([task.nEnviron,task.nPatch,nSub]);
        subRR_mean = zeros([task.nEnviron,task.nPatch,nSub]);

        for iS = 1:nSub
            for iP = 1:task.nPatch
                for iE = 1:task.nEnviron
                    subLT_mean(iE,iP,iS) = mean(trialLeaveT.leaveT(trialLeaveT.sub == iS & trialLeaveT.env == iE & trialLeaveT.patch == iP));
                    subLT_sd(iE,iP,iS) = std(trialLeaveT.leaveT(trialLeaveT.sub == iS & trialLeaveT.env == iE & trialLeaveT.patch == iP));
                    subRR_mean(iE,iP,iS) = 0;

                end
            end
        end

    case 'kane'
        data = readtable('~/Dropbox/foraging/raw_data/replication_datasets/kane/kane2019-rats-fig-1-data.csv');
        % filter to only study we need
        data = data(contains(data.Experiment,'Travel'),:);

        subID = unique(unique(data.Subject));
        nSub = numel(subID);

        subLT_mean = zeros([task.nEnviron,task.nPatch,nSub]);
        subLT_sd = zeros([task.nEnviron,task.nPatch,nSub]);
        subRR_mean = zeros([task.nEnviron,task.nPatch,nSub]);

        for iP = 1:task.nPatch
            for iE = 1:task.nEnviron
                for iS = 1:nSub
                    subjData = data(data.Subject == subID(iS), :);
                    % extract state where leave decision made, for each patch x env
                    % combination
                    subLT_mean(iE,iP,iS) = mean(subjData.StateInPatch(subjData.Decision == 1 & subjData.Travel == task.travelTime(iE)*10 & subjData.startVolume == task.r0(iP)/1000));
                    subLT_sd(iE,iP,iS) = std(subjData.StateInPatch(subjData.Decision == 1 & subjData.Travel == task.travelTime(iE)*10 & subjData.startVolume == task.r0(iP)/1000));
                    subRR_mean(iE,iP,iS) = mean(subjData.RewardVolume_mL(subjData.Decision == 1 & subjData.Travel == task.travelTime(iE)*10 & subjData.startVolume == task.r0(iP)/1000));
                end
            end
        end
end
