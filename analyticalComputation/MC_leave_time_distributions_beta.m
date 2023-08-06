% script to Monte Carlo generate leave-time distributions for different
% action selection policies - softmax and epsilon-softmax
% Sweep over range of beta parameter, for 3 different patch types (r0)
%
% Mark Humphries 13 May 2023.
% Emma Scholey updated - latest 26 July 2023

clearvars
close all;

study = 'kane'; % dataset to investigate. Options are:
% 'leheron' - human patch foraging
% 'contrerashuerta - human patch foraging
% 'kane' - rodent patch foraging

% select patch decay parameters depending on study
switch study
    case 'leheron'
        alpha = 0.075; % decay rate
        r0 = [32.5 45 57.5]; % initial yield
        t_step = 1;  % time-step at which to calculate estimates of E and VAR

    case 'contrerashuerta'
        alpha = 0.11;
        r0 = [34.5, 57.5];
        t_step = 1;

    case 'kane'
        alpha = 8;
        r0 = [60, 90, 120];
        t_step = 1;

end

action_policy = 'softmax'; % options are softmax or e-softmax

switch action_policy
    case 'softmax'
        epsilon = 0;
    case 'e-softmax'
        epsilon = 0.01;
end

number_samples = 250;  % number of Monte Carlo samples per parameter pair

beta = logspace(-3,0, 100);  % space of softmax temperatures to calculate;
% maximum of beta=1 here, as actual temperature is
% beta * r0

t_max = 100;   % maximum time in patch (for explicit calculations)
n_steps = round(t_max / t_step);  % number of time-steps

%% MC sampling
t_series = (t_step* (1:n_steps))';  % the sequence of time-steps tested

sample_leave_time = zeros(numel(beta),numel(r0),number_samples);
E_leave = zeros(numel(beta),numel(r0)); VAR_leave = zeros(numel(beta),numel(r0));
for iB = 1:numel(beta)
    for iR = 1:numel(r0)
        %% calculate expected probability of leaving on each time-step
        f_leave = zeros(n_steps,1); p = zeros(n_steps,1); r = zeros(n_steps,1);
        for n = 1:n_steps
            time_now = n * t_step;                          % what is actual time on this time-step?

            switch study
                case 'kane'
                    r(n) = reward_at_t(time_now,r0(iR),alpha,'linear');  % reward on that time-step
                otherwise
                    r(n) = reward_at_t(time_now,r0(iR),alpha, 'exp');  % reward on that time-step
            end

            p(n) = p_leave_softmax(r(n),beta(iB));          % probability on that time-step

            % probability of leaving on that time-step is the product of
            % all probabilities of *not leaving* up till the current trial and the probability of
            % leaving on this trial

            if n == 1
                f_leave(n) = p(n);
            else
                f_leave(n) = prod(1-p(1:n-1)) * p(n);
            end
        end

        E_leave(iB,iR) = sum(t_series.*f_leave);
        VAR_leave(iB,iR) = sum((t_series - E_leave(iB,iR)).^2 .* f_leave);

        %% Monte Carlo sampling of distribution
        for iSample = 1:number_samples
            blnLeave = 0; n = 1;
            while ~blnLeave
                time_now = n * t_step;                          % what is actual time on this time-step?

                switch study
                    case 'kane'
                        r = reward_at_t(time_now,r0(iR),alpha,'linear');  % reward on that time-step
                    otherwise
                        r = reward_at_t(time_now,r0(iR),alpha, 'exp');  % reward on that time-step
                end

                p = p_leave_softmax(r,beta(iB));          % probability on that time-step

                blnLeave = rand <= p;
                n = n + 1;
            end
            sample_leave_time(iB,iR,iSample) = time_now;
        end
    end
end

%% derived statistics from calculations
SD_leave = sqrt(VAR_leave);
CV_leave = E_leave ./ SD_leave;

%% statistics from MC simulations
E_MC_leave = mean(sample_leave_time,3);
SD_MC_leave = std(sample_leave_time,[],3);

%% load and summarise real data

switch study

    case 'leheron'
        load('~/Dropbox/foraging/raw_data/summary/young_variables/t_young.mat') % load real data

        nEnv = 2;
        nPatch = 3;
        nPatchDiff = 3;

        patchNames = {'low', 'medium', 'high'};
        nSub = numel(unique(t_young.subj));

        subLT_mean = zeros([nEnv,nPatch,nSub]);
        subLT_sd = zeros([nEnv,nPatch,nSub]);

        for iS = 1:nSub
            for iP = 1:nPatch
                for iE = 1:nEnv
                    subLT_mean(iE,iP,iS) = mean(t_young.leaveT(t_young.subj == iS & t_young.env == iE & t_young.patch == iP));
                    subLT_sd(iE,iP,iS) = std(t_young.leaveT(t_young.subj == iS & t_young.env == iE & t_young.patch == iP));
                end
            end
        end

        % Optimal MVT
        %ixRich = find(E_leave(:,2) > 9.5,1,'first');
        %ixPoor = find(E_leave(:,2) > 11.5,1,'first');

    case 'contrerashuerta'
        data = readtable('~/Dropbox/foraging/raw_data/replication_datasets/contrerashuerta/contrerashuerta_exp2_LT.xlsx');
        sd_data = readtable('~/Dropbox/foraging/raw_data/replication_datasets/contrerashuerta/contrerashuerta_exp2_LT_SD.xlsx');

        nEnv = 2;
        nPatch = 2;
        patchNames = {'low', 'high'};
        nPatchDiff = 1;

        nSub = height(data);

        subLT_mean = zeros([nEnv,nPatch,nSub]);
        subLT_sd = zeros([nEnv,nPatch,nSub]);

        for iS = 1:nSub
            subLT_mean(:,:,iS) = reshape(table2array(data(iS,1:4)),[2,2])';
            subLT_sd(:,:,iS) = reshape(table2array(sd_data(iS,1:4)),[2,2])';
        end

    case 'kane'
        data = readtable('~/Dropbox/foraging/raw_data/replication_datasets/kane/kane2019-rats-fig-1-data.csv');
        % filter to only study we need
        data = data(contains(data.Experiment,'Travel'),:);

        subID = unique(unique(data.Subject));
        nSub = numel(subID);

        envType = unique(data.Travel);
        patchType = unique(data.startVolume);

        nEnv = numel(envType); % 2 environments determined by travel time between patches
        nPatch = numel(patchType); % 3 patches determined by initial yield
        patchNames = {'low', 'medium', 'high'};
        nPatchDiff = 3;

        subLT_mean = zeros([nEnv,nPatch,nSub]);
        subLT_sd = zeros([nEnv,nPatch,nSub]);

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

        % plot decay function
        figure

        for iP = 1:nPatch
            plot(patchType(iP)*1000:-8:0), hold on
        end

        xlim([1,20])
        ylabel('Patch reward rate, \mu l')
        xlabel('State in patch (1 state = 10s)')
        legend({'low yield', 'medium yield', 'high yield'})

end

% group averages
mean_leave_data = mean(subLT_mean,3);
SD_leave_data = mean(subLT_sd,3);

%% plot expected vs participant leaving times
figure

c_rich = [0.7 0.3 0.3];
c_poor = [0.3 0.3 0.7];

% find predicted leaving times for this dataset:
% find index that agrees with middle (or high) -valued patch (as it occured at same
% rate in both environments)
ixRich = find(E_leave(:,2) > mean_leave_data(1,2),1,'first');
ixPoor = find(E_leave(:,2) > mean_leave_data(2,2),1,'first');

beta_rich = beta(ixRich)
beta_poor = beta(ixPoor)

% plot model predictions
plot(E_leave(ixRich,:),'.--','Color',c_rich, 'LineWidth', 2, 'MarkerSize', 15), hold on
plot(E_leave(ixPoor,:),'.--','Color',c_poor, 'LineWidth', 2, 'MarkerSize', 15)

% plot the actual experimental data
plot(mean_leave_data(1,:),'x-', 'Color',c_rich,'LineWidth', 2, 'MarkerSize', 15)
plot(mean_leave_data(2,:),'x-', 'Color',c_poor,'LineWidth', 2, 'MarkerSize', 15)

set(gca,'XLim',[0 4],'XTick',1:nPatch,'XTickLabel',patchNames, 'YLim', [0 30])
ylabel('Expected leaving time (s)')
title(sprintf('beta rich = %f, beta poor = %f', round(beta_rich,3),round(beta_poor,3)))

legend({'Model - rich', 'Model - poor', 'Subj - rich', 'Subj - poor'})

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

%% plot expected leaving times for range of beta values, compared to 'best' betas for real data

figure
semilogx(beta,E_leave); hold on
semilogx(beta,E_MC_leave); hold on
xline(beta_rich,'Color', c_rich,'Label','rich'), xline(beta_poor, 'Color', c_poor,'Label', 'poor')
xlabel('Beta: higher = exploit')
ylabel('Expected leaving time (s)')
legend({'low yield', 'medium yield', 'high yield'})

figure
semilogx(beta,SD_leave); hold on;
semilogx(beta,SD_MC_leave);
xline(beta_rich,'Color', c_rich,'Label','rich'), xline(beta_poor, 'Color', c_poor,'Label', 'poor')
xlabel('Beta: higher = exploit')
ylabel('SD of leaving time (s)')
legend({'low yield', 'medium yield', 'high yield'})

%% SD analysis
% % question: which is bigger - change in SD of patch-type between environments,
% or between patch types in same environment?

% E and SD at participant fit
E_leave_fit = [E_leave(ixRich,:); E_leave(ixPoor,:)]; %[rich; poor]
SD_leave_fit = [SD_leave(ixRich,:); SD_leave(ixPoor,:)]; %[rich; poor]

model_change_in_SD_between_environments = SD_leave_fit(2,:) - SD_leave_fit(1,:); % poor - rich
model_change_in_SD_within_rich = nonzeros(triu(SD_leave_fit(1,:) - SD_leave_fit(1,:)'))';
model_change_in_SD_within_poor = nonzeros(triu(SD_leave_fit(2,:) - SD_leave_fit(2,:)'))';

% compare to participant data

real_change_in_SD_between_environments = SD_leave_data(2,:) - SD_leave_data(1,:); % poor - rich
real_change_in_SD_within_rich = nonzeros(triu(SD_leave_data(1,:) - SD_leave_data(1,:)'))';
real_change_in_SD_within_poor = nonzeros(triu(SD_leave_data(2,:) - SD_leave_data(2,:)'))';

% plot model vs participant SDs
figure

plot(SD_leave_fit(1,:),'.--','Color',c_rich, 'LineWidth', 2, 'MarkerSize', 15); hold on
plot(SD_leave_fit(2,:),'.--','Color',c_poor, 'LineWidth', 2,'MarkerSize', 15);
plot(SD_leave_data(1,:),'x-','Color',c_rich, 'LineWidth', 2,'MarkerSize', 15)
plot(SD_leave_data(2,:),'x-','Color',c_poor, 'LineWidth', 2,'MarkerSize', 15)
scatter(1:nPatch, squeeze(subLT_sd(1,:,:)),[],c_rich, 'filled', 'MarkerFaceAlpha', 0.2)
scatter(1:nPatch, squeeze(subLT_sd(2,:,:)),[],c_poor, 'filled', 'MarkerFaceAlpha', 0.2)

set(gca,'XLim',[0 4],'XTick',1:nPatch,'XTickLabel',patchNames)
ylabel('SD of leaving time (s)')
legend({'Model - rich', 'Model - poor', 'Data - rich', 'Data - poor'})


% plot between environment differences
figure

plot(model_change_in_SD_between_environments,'.--', 'Color', 'black', 'LineWidth', 2, 'MarkerSize', 15); hold on
plot(real_change_in_SD_between_environments,'x-', 'Color', 'black', 'LineWidth', 2,'MarkerSize', 15)
scatter(1:nPatch, squeeze(subLT_sd(2,:,:)-subLT_sd(1,:,:)),'k', 'filled', 'MarkerFaceAlpha', 0.2)

set(gca,'XLim',[0 4],'XTick',1:nPatch,'XTickLabel',patchNames)
ylabel('Change in SD leaving time (s) (poor - rich)')
legend({'Model' 'Data'})

% plot within environment differences
figure

points_real_change_in_SD_within = zeros([nSub,nPatchDiff,nEnv]);
for iE = 1:nEnv
    switch study
        case 'contrerashuerta'
            points_real_change_in_SD_within(:,1,iE) = squeeze(subLT_sd(iE,2,:)-subLT_sd(iE,1,:)); % high-low
            patchDiffNames = {'High-low'};
        otherwise
            points_real_change_in_SD_within(:,1,iE) = squeeze(subLT_sd(iE,2,:)-subLT_sd(iE,1,:)); % mid-low
            points_real_change_in_SD_within(:,2,iE) = squeeze(subLT_sd(iE,3,:)-subLT_sd(iE,1,:)); % high-low
            points_real_change_in_SD_within(:,3,iE) = squeeze(subLT_sd(iE,3,:)-subLT_sd(iE,2,:)); % high-mid
            patchDiffNames = {'Mid-Low','High-low','High-mid'};
    end
end

plot(model_change_in_SD_within_rich,'.--','Color',c_rich, 'LineWidth', 2, 'MarkerSize', 15); hold on
plot(model_change_in_SD_within_poor,'.--','Color',c_poor, 'LineWidth', 2,'MarkerSize', 15);
plot(real_change_in_SD_within_rich,'x-','Color',c_rich, 'LineWidth', 2,'MarkerSize', 15) % rich
plot(real_change_in_SD_within_poor,'x-','Color',c_poor, 'LineWidth', 2,'MarkerSize', 15) % poor
scatter(1:nPatchDiff, points_real_change_in_SD_within(:,:,1),[],c_rich, 'filled', 'MarkerFaceAlpha', 0.2)
scatter(1:nPatchDiff, points_real_change_in_SD_within(:,:,2),[],c_poor, 'filled', 'MarkerFaceAlpha', 0.2)

set(gca,'XLim',[0 4],'XTick',1:nPatchDiff,'XTickLabel',patchDiffNames)
ylabel('Change in SD leaving time (s) (patch yield)')
legend({'Model - rich', 'Model - poor', 'Data - rich', 'Data - poor'})
