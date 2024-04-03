%% load participant data and summarise into leaving times
function [subLT_mean, subLT_sd, subRR_mean] = loadLeaveT(study,task)

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

