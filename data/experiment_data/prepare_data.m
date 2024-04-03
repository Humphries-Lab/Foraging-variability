%% le heron dataset 

clear
lt = readtable('leheron_trialbytrial.csv','ReadVariableNames',true);

for iS = unique(lt.sub)'

    their_data = lt(lt.sub==iS,:);

    switchIndex = diff(their_data.env) ~= 0;
    blockSwitchIndex{iS} = [1; switchIndex]; % append with 1 (first patch in first block);
end

save('leheron_blockSwitchIndex.mat','blockSwitchIndex')

%% contreras huerta dataset
clear
lt = readtable('contrerashuerta_trialbytrial.csv','ReadVariableNames',true);

data = table;
for iS = [450:452,454:450+lt.sub(end)]
    tmp = load(sprintf('seb_raw_data/%d.mat',iS));
    subj_data=struct2table(tmp.result.data);
    subj_data = subj_data(subj_data.tEnd>0.5,:);

    subj_data = subj_data(~isoutlier(subj_data.tEnd,'mean'),:); %% don't use interquartlie range ('quartiles') as it's too strict. 
    subj_data = subj_data(subj_data.condN == 1,:); % restrict to self beneficiary only
    data = table(subj_data.block, subj_data.blVal);

    subjOrder = unique(data, 'rows','stable');
    blockOrder(iS-449,:) = subjOrder.Var2;
    switchIndex = diff(data.Var1) ~= 0;  
    blockSwitchIndex{iS-449} = [1; switchIndex]; % append with 1 (first patch in first block); 

end

blockOrder = blockOrder([1:3,5:end],:);
blockSwitchIndex = blockSwitchIndex([1:3,5:end]);

writematrix(blockOrder, 'contrerashuerta_blockOrder.csv')
save('contrerashuerta_blockSwitchIndex.mat','blockSwitchIndex')

%% kane dataset 
% get _trialbytrial in same format as leheron and CH
clear
T = readtable('kane2019-rats-fig-1-data.csv','ReadVariableNames',true);
T = T(contains(T.Experiment,'Travel'),:); % only looking at travel time experiment
T_leave = T(T.Decision == 1,:); % restrict to decision states (i.e. per patch) 

endTrialsIndex = find(T.Trial==1) - 1; % find the last trial of each session
endTrialsIndex(end+1) = numel(T.Trial); % add on last trial
T_end = T(endTrialsIndex(2:end), :);

data = table(T_leave.Subject, T_leave.startVolume, T_leave.Travel, T_leave.StateInPatch);
data.Properties.VariableNames = {'sub', 'patch', 'env', 'leaveT'};
data.leaveT = round(data.leaveT);

data.patch(data.patch == 0.06) = 1; % convert patches to [1 2 3]
data.patch(data.patch == 0.09) = 2;
data.patch(data.patch == 0.12) = 3;

data.env(data.env == 10) = 1;
data.env(data.env == 30) = 2;

% switch index and block order
subjectNum = unique(T_leave.Subject);

for iS = 1:numel(subjectNum)
    tmp = T_leave.Date(T_leave.Subject == subjectNum(iS));
    switchIndex = diff(tmp) ~= 0; 
    blockSwitchIndex{iS} = [1; switchIndex]; % append with 1 (first patch in first block); 

    tmp = data.env(data.sub == subjectNum(iS));
    tmp = tmp(logical(blockSwitchIndex{iS}));
    blockOrder(iS,:) = tmp';
end
save('kane_blockSwitchIndex.mat','blockSwitchIndex')
writematrix(blockOrder, 'kane_blockOrder.csv')

% calculate experienced avgRR for each rodent, and optimal rates
tmp = unique([T_leave.Subject, T_leave.CumulativeRate, T_leave.Travel], 'rows', 'stable');
travelTime = [10 30];
for iE = 1:2
    optAvgRR(iE) = mean(tmp(tmp(:,3)==travelTime(iE),2)); % MVT optimal rates (I think) 
end

for iE = 1:2
    for iS = 1:numel(subjectNum)
        % find experienced avgRR for each rodent
        tmp = T_end(T_end.Subject == subjectNum(iS),:);
        averageRR = tmp.CumulativeReward_mL./tmp.CumulativeTime;
        experiencedAvgRR(iS,iE) = mean(averageRR(tmp.TravelTime == travelTime(iE)));
    end

end
save('kane_MVT_optimal_rates.mat', 'optAvgRR')
writematrix(experiencedAvgRR, 'kane_experiencedAvgRR.csv')


for iS = 1:numel(subjectNum)
    data.sub(data.sub==subjectNum(iS)) = iS;
end
writetable(data, 'kane_trialbytrial.csv')

%% barack

clear
load('data_adhd.mat')

travelTimes = [1 5];

% exclude subj 370 - only did one environment
d(370) = [];

nSubj = length(d);
nEnv = numel(travelTimes);

for iS = 1:nSubj
    for iE = 1:nEnv
        env_leaves = d(iS).leaves(d(iS).tt == travelTimes(iE));
        find_leave_times = find(env_leaves);
        leave_times = [find_leave_times(1);diff(find_leave_times)]-1;
        tmp{iE}.leaveT = leave_times';
        tmp{iE}.patch = repelem(1,numel(tmp{iE}.leaveT));
        tmp{iE}.env = repelem(iE,numel(tmp{iE}.leaveT)); 
        % the -1 is because the 'leave' action is on that time point i.e. no reward immediately
        tmp{iE}.sub = repelem(iS,numel(tmp{iE}.leaveT));
    end
        table(tmp{1}.leaveT, tmp{1}.patch, tmp{1}.env, tmp{1}.sub)

    [struct2table(tmp{1}), struct2table(tmp(2))]
            patchData{iS} = struct2table([tmp{1}, tmp{2}])

end


lt = readtable('barack_trialbytrial.csv','ReadVariableNames',true);

for iS = unique(lt.sub)'

    their_data = lt(lt.sub==iS,:);

    switchIndex = diff(their_data.env) ~= 0;
    blockSwitchIndex{iS} = [1; switchIndex]; % append with 1 (first patch in first block);
end

save('leheron_blockSwitchIndex.mat','blockSwitchIndex')
