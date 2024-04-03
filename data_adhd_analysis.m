%% script to check SD of David Barack data with linear decay 

clear
load('../raw_data/data_adhd.mat')

travelTimes = [1 5];

% exclude subj 370 - only did one environment
d(370) = [];

nSubj = length(d);
nEnv = numel(travelTimes);

mean_leaveT = zeros(nSubj, nEnv); std_leaveT = zeros(nSubj, nEnv);

for iS = 1:nSubj
    for iE = 1:nEnv
        env_leaves = d(iS).leaves(d(iS).tt == travelTimes(iE));
        find_leave_times = find(env_leaves);
        leave_times = [find_leave_times(1);diff(find_leave_times)]-1; % the -1 is because the 'leave' action is on that time point i.e. no reward immediately
        mean_leaveT(iS,iE) = mean(leave_times);
        std_leaveT(iS,iE) = std(leave_times);
        nPatches_visited(iS,iE) = sum(d(iS).leaves(d(iS).tt == travelTimes(iE)));
    end
end

%% quick plots 
close all

tmp_set1 = brewermap(9,'Set1');
color.rich = tmp_set1(2,:); % green
color.poor = tmp_set1(3,:); % blue


% exclude people with 0 variability - just counting 
index1 = find(std_leaveT(:,1)==0);std_leaveT(index1,:) = []; mean_leaveT(index1,:) = []; 
index2 = find(std_leaveT(:,2)==0);std_leaveT(index2,:) = []; mean_leaveT(index2,:) = [];

% summarise data
subjectMean = median(mean_leaveT); subjectMean_SEM = std(mean_leaveT)./nSubj;
subjectSD = median(std_leaveT); subjectSD_SEM = std(std_leaveT)./nSubj;

% mean panel
figure

% plot the actual experimental data
errorbar(1, subjectMean(1),subjectMean_SEM(1), 'Color', [0 0 0], 'LineWidth',3), hold on
errorbar(2, subjectMean(2),subjectMean_SEM(2), 'Color', [0 0 0], 'LineWidth',3)
set(gca,'XLim',[0, 3],'XTick',1:2,'XTickLabel',{'rich', 'poor'})
ylabel('Median stay decisions')
set(findall(gcf,'-property','FontSize'),'FontSize',18)

% SD panel
figure

% plot the actual experimental data
errorbar(1, subjectSD(1),subjectSD_SEM(1), 'Color', [0 0 0], 'LineWidth',3), hold on 
errorbar(2, subjectSD(2),subjectSD_SEM(2), 'Color', [0 0 0], 'LineWidth',3)
set(gca,'XLim',[0, 3],'XTick',1:2,'XTickLabel',{'rich', 'poor'})
ylabel('Median SD stay decisions')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
