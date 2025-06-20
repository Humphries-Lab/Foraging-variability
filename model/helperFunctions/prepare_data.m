%% Data Preparation for Foraging Variability project
% code to prepare raw data into derived data files for model fitting and
% analysis

%% le heron dataset 
% Emma Scholey 16 January 2015

cd('../../data/experiment_data/')

clear
lt = readtable('./leheron/leheron_trialbytrial.csv','ReadVariableNames',true);

for iS = unique(lt.sub)'

    their_data = lt(lt.sub==iS,:);

    switchIndex = diff(their_data.env) ~= 0;
    blockSwitchIndex{iS} = [1; switchIndex]; % append with 1 (first patch in first block);
end

blockOrder = readmatrix("leheron/leheron_blockOrder.csv");
experiencedAvgRR = readmatrix("leheron/leheron_experiencedAvgRR.csv");

subj_var.blockOrder = blockOrder;
subj_var.blockSwitchIndex = blockSwitchIndex;
subj_var.experiencedAvgRR = experiencedAvgRR; 
save('leheron/leheron_subj_var.mat','subj_var')

%% contreras huerta dataset
clear
lt = readtable('contrerashuerta/contrerashuerta_trialbytrial.csv','ReadVariableNames',true);

data = table;

subjNumber = [450:452,454:450+lt.sub(end)];

for iS = 1:length(subjNumber)
    lt_subj_data = lt(lt.sub == iS & lt.ben == 1,:);

    tmp = load(sprintf('../../../raw_data/contrerashuerta/contrerashuerta_raw_data/%d.mat',subjNumber(iS)));
    subj_data=struct2table(tmp.result.data);
    subj_data = subj_data(subj_data.condN == 1,:); % restrict to self beneficiary only


    % filter the trials from raw data to exclude the ones in lt by Seb
    % (e.g. outlier or lower than 0.5). 

    subj_data = subj_data(ismember(subj_data.tEnd, lt_subj_data.leaveT),:);

    data = table(subj_data.block, subj_data.blVal);

    subjOrder = unique(data, 'rows','stable');
    blockOrder(iS,:) = subjOrder.Var2;
    switchIndex = diff(data.Var1) ~= 0;  
    blockSwitchIndex{iS} = [1; switchIndex]; % append with 1 (first patch in first block); 

end

experiencedAvgRR = readmatrix("contrerashuerta/contrerashuerta_experiencedAvgRR.csv");

subj_var.blockOrder = blockOrder;
subj_var.blockSwitchIndex = blockSwitchIndex;
subj_var.experiencedAvgRR = experiencedAvgRR; 

save('contrerashuerta/contrerashuerta_subj_var.mat','subj_var')

%% kane dataset 
% get _trialbytrial in same format as leheron and CH
clear
T = readtable('kane/kane2019-rats-fig-1-data.csv','ReadVariableNames',true);
T = T(contains(T.Experiment,'Travel'),:); % only looking at travel time experiment

T = sortrows(T,[2,4]); % sort by subject (column 2) and date (column 4)

T_leave = T(T.Decision == 1,:); % restrict to decision states (i.e. per patch) 

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

travelTime = [10 30];

for iS = 1:numel(subjectNum)
    tmp = T_leave.Date(T_leave.Subject == subjectNum(iS));
    switchIndex = diff(tmp) ~= 0; 
    blockSwitchIndex{iS} = [1; switchIndex]; % append with 1 (first patch in first block); 

    tmp = data.env(data.sub == subjectNum(iS));
    tmp = tmp(logical(blockSwitchIndex{iS}));
    blockOrder(iS,:) = tmp';

    % calculate experienced avgRR for each rat
    tmp = T_leave(T_leave.Subject == subjectNum(iS),:);
    tmp = tmp(logical([switchIndex;1]),:); % append last trial on switch index, to account for last trial
    averageRR = tmp.CumulativeReward_mL./tmp.CumulativeTime;

    for iE = 1:2
        experiencedAvgRR(iS,iE) = 10000 * mean(averageRR(tmp.TravelTime == travelTime(iE)));
    end
end
 
% calculate optimal rates for each rat
tmp = unique([T_leave.Subject, T_leave.CumulativeRate, T_leave.Travel], 'rows', 'stable');
for iE = 1:2
    optAvgRR(iE) = mean(tmp(tmp(:,3)==travelTime(iE),2))*10000; % MVT optimal rates
end

subj_var.blockOrder = blockOrder;
subj_var.blockSwitchIndex = blockSwitchIndex;
subj_var.experiencedAvgRR = experiencedAvgRR; 
subj_var.optAvgRR = optAvgRR;

save('kane/kane_subj_var.mat', 'subj_var')

% make subject numbers 1:8 rather than pID
for iS = 1:numel(subjectNum)
    data.sub(data.sub==subjectNum(iS)) = iS;
end
 %writetable(data, 'kane_trialbytrial.csv')

