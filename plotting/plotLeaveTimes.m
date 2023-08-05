%% Plot participant data
function plotLeaveTimes(study,task,modelLT)

switch study

    case 'leheron'
        load('~/Dropbox/foraging/raw_data/summary/young_variables/t_young.mat') % load real data

        nSub = numel(unique(t_young.subj));

        subLT_mean = zeros([task.nEnviron,task.nPatch,nSub]);
        subLT_sd = zeros([task.nEnviron,task.nPatch,nSub]);

        for iS = 1:nSub
            for iP = 1:task.nPatch
                for iE = 1:task.nEnviron
                    subLT_mean(iE,iP,iS) = mean(t_young.leaveT(t_young.subj == iS & t_young.env == iE & t_young.patch == iP));
                    subLT_sd(iE,iP,iS) = std(t_young.leaveT(t_young.subj == iS & t_young.env == iE & t_young.patch == iP));
                end
            end
        end

    case 'contrerashuerta'
        data = readtable('~/Dropbox/foraging/raw_data/replication_datasets/contrerashuerta_exp2_LT.xlsx');
        sd_data = readtable('~/Dropbox/foraging/raw_data/replication_datasets/contrerashuerta_exp2_LT_SD.xlsx');

        nSub = height(data);

        subLT_mean = zeros([task.nEnviron,task.nPatch,nSub]);
        subLT_sd = zeros([task.nEnviron,task.nPatch,nSub]);

        for iS = 1:nSub
            subLT_mean(:,:,iS) = reshape(table2array(data(iS,1:4)),[2,2])';
            subLT_sd(:,:,iS) = reshape(table2array(sd_data(iS,1:4)),[2,2])';
        end

    case 'kane'
        data = readtable('~/Dropbox/foraging/raw_data/replication_datasets/kane2019-rats-fig-1-data.csv');
        % filter to only study we need
        data = data(contains(data.Experiment,'Travel'),:);

        subID = unique(unique(data.Subject));
        nSub = numel(subID);

        subLT_mean = zeros([task.nEnviron,task.nPatch,nSub]);
        subLT_sd = zeros([task.nEnviron,task.nPatch,nSub]);

        for iP = 1:nPatch
            for iE = 1:nEnv
                for iS = 1:nSub
                    subjData = data(data.Subject == subID(iS), :);

                    % extract state where leave decision made, for each patch x env
                    % combination
                    subLT_mean(iE,iP,iS) = mean(subjData.StateInPatch(subjData.Decision == 1 & subjData.Travel == envType(iE) & subjData.startVolume == patchType(iP)));
                    subLT_sd(iE,iP,iS) = std(subjData.StateInPatch(subjData.Decision == 1 & subjData.Travel == envType(iE) & subjData.startVolume == patchType(iP)));
                end
            end
        end
end

% group averages
subj_leave_data = mean(subLT_mean,3);

model_leave_data = mean(modelLT,3);

%% plot participant leaving times
figure

c_rich = [0.7 0.3 0.3];
c_poor = [0.3 0.3 0.7];

% plot the simulated data
plot(model_leave_data(1,:),'.--', 'Color',c_rich,'LineWidth', 2, 'MarkerSize', 15), hold on
plot(model_leave_data(2,:),'.--', 'Color',c_poor,'LineWidth', 2, 'MarkerSize', 15)

% plot the actual experimental data
plot(subj_leave_data(1,:),'.-', 'Color',c_rich,'LineWidth', 2, 'MarkerSize', 15)
plot(subj_leave_data(2,:),'.-', 'Color',c_poor,'LineWidth', 2, 'MarkerSize', 15)

set(gca,'XLim',[0 4],'XTick',1:task.nPatch,'XTickLabel',task.patchNames, 'YLim', [0 30])
ylabel('Leaving time (s)')

legend({'model - rich', 'model - poor', 'data - rich', 'data - poor'})

% overlay optimal MVT predictions
% % overlay MVT optimal times
% RR_Leave = [21.8678 18.5632]; % this is the average background RR at which subjects should leave, for Rich and Poor env... -
% ... respectively.
%
% clear optLT
% A=[32.5 45 57.5];a=0.075;
% for e=1:2 %each env
%     for p=1:3 % each patch type
%         optLT(e,p)=(log(RR_Leave(e)/A(p)))/-a;
%     end
% end
% plot(optLT(1,:),'--','Color',c_rich)
% plot(optLT(2,:),'--','Color',c_poor)