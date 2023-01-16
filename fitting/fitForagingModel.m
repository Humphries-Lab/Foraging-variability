% fitting RL model functions to foraging data
% Emma Scholey 9 Jun 2022
clear
close all;
addpath(genpath('~/Dropbox/foraging/code'))

model = 25;

SetUpEnviron % get the environment parameters

% fmincon options
lowerBounds = [0,0,0];  % [alpha_Q, alpha_rho, beta] % parameter bounds
upperBounds = [1,1,100];  % arbitrary upper bound on beta to stop pathological behaviour
options = optimoptions('fmincon','Display','none'); % don't display

% get their data
%subj = [1:23 25:40];
subj = 1:4;
load ~/Dropbox/foraging/raw_data/summary/young_variables/t_young.mat
group = {'Young_HC_ReDo/yc%g_forage.mat', subj};

for k = 1:length(subj) % for each subject
    for b = 1:2 % for each block [rich poor]
    SubjLT = t_young(t_young.subj == k,:); % extract their summarised leaving times
    Data{k}{b} = PrepSubjData(SubjLT, b, Env); % transform into state actions ready for fitting
    end
end

clear group k b SubjLT t_young 

% initialise
numSubjects = size(subj,2);
numParams = size(lowerBounds,2);

%% grid search for subjects - marginal likelilood distributions
% range of values for grid search of starting points for search
AlphaVector = 0:0.05:1;
BetaVector = 0:2.5:50;

for k = 1:numSubjects
    k
    for block = 1:2 % for each environment: [rich, poor]
        [MinNLLEval{k}{block}, ParamsMinNLLEval{k}{block}, NLLEval{k}{block}] = GridSearchForaging(Data{k}{block}, Env, AlphaVector, BetaVector); % change this depending on model to test
    end
end

% plot marginal likelihoods across grid for each subject
plotMarginalLikelihoods(NLLEval, AlphaVector,BetaVector)


%% fitting for each person in group with different starting points
nStarts = 50; % how many starting points to avoid local minima

each_params_fitted = zeros(numSubjects,numParams,nStarts);
each_NLL_eval = zeros(numSubjects,nStarts);

for k = 1:numSubjects
    for block = 1:2
        for ii = 1:nStarts
            [SubjMinNLLFitParams, SubjMinNLLEval, Hessian, StdErrorMinNLLFitParams] = Fit_M25_RLOn(Data(k),Env);
        end
    end
end

for ii = 1:nStarts
    params0(ii,:) = [rand, rand, exprnd(2)] % randomly select parameters each time
    for k = 1:numSubjects
            f = @(x)NLL_subject_foraging_RL_on(x,SubjData{k},Patch_Order{k}, environ,model_table, model);
        [each_params_fitted(k,:,ii),each_NLL_eval(k,ii), ~, ~, ~, ~, hessian(:, :, k, ii)] = ...
            fmincon(f,params0(ii,:),[],[],[],[],lowerBounds,upperBounds, [], options);

        % calculate standard error on the parameter fits using the hessian
        % matrix
        tmp = sqrt(diag(inv(hessian(:,:,k,ii))))'; % The diagonal terms of H−1 correspond to variances for each parameter separately, and their square roots measure one standard error on the parameter. - Daw 2011
        tmp(~isreal(tmp),:) = nan;
        SE_params_fitted(k,:,ii) = tmp;
    end
end

%% Each subject's best fitting parameters

% find the parameters corresponding to the NLL minima
minNLL = zeros(numSubjects,1);
minNLL_params_fitted = zeros(numSubjects,numParams);
minNLL_SE_params_fitted = zeros(numSubjects,numParams);
BIC = zeros(numSubjects,1);
AIC = zeros(numSubjects,1);

% get the number of data points for each subject (value of params doesn't
% matter
[~, numObservations] = NLL_group_foraging([0 0 1],SubjData, Patch_Order, environ, model_table, model);

for k = 1:numSubjects
    % for each subject, find their minimum NLL across grid of starting points
    their_NLL = squeeze(each_NLL_eval(k,:));
    minNLL(k) = min(their_NLL);   % minimum negative log likelihood over all starting positions
    ix = find(their_NLL == minNLL(k));    % indices of location of minimum

    % get the corresponding fitted parameter values. Note: if more than one starting location converges on same parameters then will have same NLL, so use first
    minNLL_params_fitted(k,:) = squeeze(each_params_fitted(k,:,ix));

    % get the corresponding SE for the fitted parameter.
    minNLL_SE_params_fitted(k,:) = squeeze(SE_params_fitted(k,:,ix));

    % BIC for each person
    BIC(k) = numParams * log(numObservations(k)) + 2*minNLL(k);
    AIC(k) = 2/numObservations(k) * minNLL(k) + 2 * numParams/numObservations(k);
end

names = {'rich', 'poor'};
save_name = sprintf('~/Dropbox/foraging-project/results/M%d/fitting_results_%s', model, names{block});
save(save_name, 'AIC', 'BIC', 'minNLL', 'each_params_fitted', 'minNLL_params_fitted', 'minNLL_SE_params_fitted')

%% plots
% parameter distribution for each subject, based on different starting
% points for fmincon
close all
figure;
tl = tiledlayout('flow', 'TileSpacing', 'Compact');
for k = 1:4
    ax = nexttile;
    histogram(squeeze(each_params_fitted(k,1,:))); % alpha patch
    title(sprintf('Subject %d', k))
    xlabel('alpha patch - fitted', 'FontSize', 12)
end
title(tl, 'distribution of fit parameter values - alpha patch', 'FontSize', 16);

figure;
tl = tiledlayout('flow', 'TileSpacing', 'Compact');
for k = 1:4
    nexttile;
    histogram(each_params_fitted(k,2,:), 'FaceColor',[0.8500 0.3250 0.0980]);  % alpha rho
    title(sprintf('Subject %d', k))
    xlabel('alpha rho - fitted', 'FontSize', 12)
end
title(tl, 'distribution of fit parameter values - alpha rho', 'FontSize', 16);

figure;
tl = tiledlayout('flow', 'TileSpacing', 'Compact');
for k = 1:4
    nexttile;
    histogram(each_params_fitted(k,3,:)); % beta
    title(sprintf('Subject %d', k))
    xlabel('beta - fitted', 'FontSize', 12)
end
title(tl, 'distribution of fit parameter values - beta', 'FontSize', 16);

% plot the best parameters for each person
m = median(minNLL_params_fitted);

figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');

nexttile
plot(minNLL_params_fitted(:,1),minNLL_params_fitted(:,2),'o','MarkerSize',10); hold on
plot(m(1),m(2),'+','MarkerSize',15)
xlabel('Q/patch learning rate (\alpha)','FontSize', 12);
ylabel('\rho learning rate (\alpha)','FontSize', 12);

nexttile
plot(minNLL_params_fitted(:,1),minNLL_params_fitted(:,3),'o','MarkerSize',10); hold on
plot(m(1),m(3),'+','MarkerSize',15)
xlabel('Q/patch learning rate (\alpha)','FontSize', 12);
ylabel('Inverse temperature (\beta)','FontSize', 12);

nexttile
plot(minNLL_params_fitted(:,2),minNLL_params_fitted(:,3),'o','MarkerSize',10); hold on
plot(m(2),m(3),'+','MarkerSize',15)
xlabel('\rho learning rate (\alpha)', 'FontSize', 12);
ylabel('Inverse temperature (\beta)','FontSize', 12);
title(tl, 'Distribution of best fitting parameters', 'FontSize', 16);

leg = legend({'Subject fit', 'Group Median'}, 'FontSize', 14);
leg.Layout.Tile = 'north';

nexttile
histogram(minNLL)
title('Distribution of subject min NLL','FontSize', 12)

% inspect standard errors of the best fit parameters
figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');
nexttile
errorbar(subj, minNLL_params_fitted(:,1), minNLL_SE_params_fitted(:,1), 'LineStyle', 'none','Marker', 'o')
title('Q/patch learning rate (\alpha)','FontSize', 12);
xlabel('Subject number', 'FontSize', 12)
ylabel('Standard error', 'FontSize', 12)
nexttile
errorbar(subj, minNLL_params_fitted(:,2), minNLL_SE_params_fitted(:,2), 'LineStyle', 'none','Marker', 'o')
title('rho learning rate (\alpha)','FontSize', 12);
xlabel('Subject number', 'FontSize', 12)
ylabel('Standard error', 'FontSize', 12)
nexttile
errorbar(subj, minNLL_params_fitted(:,3), minNLL_SE_params_fitted(:,3), 'LineStyle', 'none','Marker', 'o')
title('beta (\beta)', 'FontSize', 12);
xlabel('Subject number', 'FontSize', 12)
ylabel('Standard error', 'FontSize', 12)
title(tl, 'Distribution of standard errors of best fit parameters', 'FontSize', 16);

% plot to show standard errors, excluding those with massive SE
figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');
nexttile
iSE = minNLL_SE_params_fitted(:,1)<1;
errorbar(subj(iSE), minNLL_params_fitted(iSE,1), minNLL_SE_params_fitted(iSE,1), 'LineStyle', 'none', 'Marker', 'o')
title('Q/patch learning rate (\alpha)','FontSize', 12);
xlabel('Subject number', 'FontSize', 12)
ylabel('Standard error', 'FontSize', 12)

nexttile
iSE = minNLL_SE_params_fitted(:,2)<1;
errorbar(subj(iSE), minNLL_params_fitted(iSE,2), minNLL_SE_params_fitted(iSE,2), 'LineStyle', 'none','Marker', 'o')
title('rho learning rate (\alpha)','FontSize', 12);
xlabel('Subject number', 'FontSize', 12)
ylabel('Standard error', 'FontSize', 12)
nexttile
iSE = minNLL_SE_params_fitted(:,3)<20;
errorbar(subj(iSE), minNLL_params_fitted(iSE,3), minNLL_SE_params_fitted(iSE,3), 'LineStyle', 'none','Marker', 'o')
title('beta (\beta)', 'FontSize', 12);
xlabel('Subject number', 'FontSize', 12)
ylabel('Standard error', 'FontSize', 12)
title(tl, 'Distribution of standard errors of best fit parameters (excluding outliers)', 'FontSize', 16);

% plot distance from start to fit parameters

figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');

for k = 1:4 % just test on a few participants
    nexttile;
    scatter3(squeeze(each_params_fitted(k,1,:)),squeeze(each_params_fitted(k,2,:)),squeeze(each_params_fitted(k,3,:)), 50, squeeze(each_params_fitted(k,3,:)), 'filled');
    hold on;
    scatter3(params0(:,1),params0(:,2),params0(:,3), 50, params0(:,3), "r", 'filled');
    colormap("winter");
    colorbar;
    title(sprintf('Subject %d', k))
    xlabel('alpha patch', 'FontSize', 12)
    ylabel('alpha rho', 'FontSize', 12)
    zlabel('beta', 'FontSize', 12)
    % plot lines between the start and fit points
    for ii = 1:nStarts
        v1 = [params0(ii,1),params0(ii,2),params0(ii,3)];
        v2 = [squeeze(each_params_fitted(k,1,ii)),squeeze(each_params_fitted(k,2,ii)),squeeze(each_params_fitted(k,3,ii))];
        v = [v2;v1];
        plot3(v(:,1), v(:,2),v(:,3),'r', 'LineStyle', '--')
    end
end
title(tl, 'Distance from start to fit parameters', 'FontSize', 16)
leg = legend({'Fit parameters', 'Start parameters'}, 'FontSize', 14);
leg.Layout.Tile = 'north';


%% parameter recovery for range of values

Block = 1;

nSim = 250; % number of iterations

SimParams = zeros([numParams,nSim]);
FitParams = zeros([numParams,nSim]);

for ii = 1:nSim
    %TestParams = [0.001+rand(1,1)*(0.6-0.001), 0.001+rand(1,1)*(0.15-0.001), exprnd(10)]; % randomly select parameters each time
    TestParams = [rand rand exprnd(10)]
    % simulate data
    out = simulateM25_RLOn(TestParams, Block, Env); % change this depending on model to test
    SimData.Action = out.Action;
    SimData.PatchOrder = Env.PatchOrder{Block};
    SimData.NextPatch = out.NextPatch; 
    % recover data
    f = @(x)NLL_M25_RLOn(x, Env, SimData);

    for nStarts = 1:10 % do a few starting points to try and avoid global minima
        %params0 = [0.001+rand(1,1)*(0.6-0.001), 0.001+rand(1,1)*(0.15-0.001), exprnd(10)]; % starting parameters
        params0 = [rand rand exprnd(10)];
        [SimParamsFitted(nStarts,:),NLLEval(nStarts)] = fmincon(f,params0,[],[],[],[],lowerBounds,upperBounds, [], options); % start search from random location each time
    end

    % get the global minima
    ix = find(NLLEval == min(NLLEval));    % indices of location of minimum negative log likelihood over all starting positions
    % get the corresponding fitted parameter values. Note: if more than one starting location converges on same parameters then will have same NLL, so use first
    minNLLSimParamsFitted = SimParamsFitted(ix,:);

    SimParams(1,ii) = TestParams(1); % alpha patch rr simulated
    SimParams(2,ii) = TestParams(2); % alpha rho simualted
    SimParams(3,ii) = TestParams(3); % beta simulated

    FitParams(1,ii) = minNLLSimParamsFitted(1); % alpha patch rr fit
    FitParams(2,ii) = minNLLSimParamsFitted(2); % alpha rho fit
    FitParams(3,ii) = minNLLSimParamsFitted(3); % beta fit
% 
%     % get their Q_stay table (based on the fit parameter values)
%     if strcmp(model_table.policy{model}, 'on')
%         [~,~,Q_stay_tables{ii},~,~,~] = NLL_subject_foraging_RL_on(minNLLSimParamsFitted,sim_foraging,patch_order,environ,model_table,model);
%     elseif strcmp(model_table.policy{model}, 'off')
%         [~,~,Q_stay_tables{ii},~,~,~] = NLL_subject_foraging_RL_off(minNLLSimParamsFitted,sim_foraging,patch_order,environ,model_table,model);
%     end

end

%% plot correlation - simulated value versus actual value
close all
names = {'learning rate - Q/patch' 'learning rate - rho' 'softmax temperature'};
symbols = {'\alpha' '\alpha' '\beta'};
figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');

% for each parameter, plot sim vs fit 
for i= 1:size(SimParams,1)
    ax = nexttile;
    plot(SimParams(i,:), FitParams(i,:), 'o', 'markersize', 8, 'linewidth', 1)
    corrcoef(SimParams(i,:),FitParams(i,:))
    %corr(SimParams(i,:)',FitParams(i,:)', 'Type', 'Spearman')
    title(sprintf('%s', names{i}))
    xlabel(sprintf('simulated %s', symbols{i}))
    ylabel(sprintf('fit %s', symbols{i}))
end
title(tl, 'Parameter recovery - correlations between simulated and fit parameters')
set(ax,'xscale', 'log', 'yscale' ,'log')

corrcoef(log(SimParams(3,:)),log(FitParams(3,:)))

figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');
nexttile;
plot(FitParams(1,:), FitParams(2,:), 'o', 'markersize', 8, 'linewidth', 1)
title('fit alpha Q/patch vs fit alpha rho')
xlabel('fit alpha Q/patch')
ylabel('fit alpha rho')

ax = nexttile;
plot(FitParams(2,:), FitParams(3,:), 'o', 'markersize', 8, 'linewidth', 1)
title('fit alpha rho vs fit beta')
xlabel('fit alpha rho')
ylabel('fit beta')
set(ax,'yscale' ,'log')

ax = nexttile;
plot(FitParams(1,:), FitParams(3,:), 'o', 'markersize', 8, 'linewidth', 1)
title('fit alpha Q/patch vs fit beta')
xlabel('fit alpha Q/patch')
ylabel('fit beta')
set(ax,'yscale' ,'log')
title(tl, 'Parameter trade-offs - correlations between fit parameters')

%% Checks for parameter recovery - extreme states 
for ii = 1:nSim
    max_Q_length(ii) = max(sum(Q_stay{ii} ~= 0));
end

residuals_alpha_patch = sqrt((FitParams(1,:) - SimParams(1,:)).^2);
residuals_alpha_rho = sqrt((FitParams(2,:) - SimParams(2,:)).^2);
residuals_beta = sqrt((FitParams(3,:) - SimParams(3,:)).^2);

% lower distance should correlate with shorter Q length
corrcoef(log(max_Q_length), log(residuals_alpha_patch))
corr(max_Q_length', residuals_alpha_patch', 'Type', 'Spearman')

corrcoef(log(max_Q_length), log(residuals_alpha_rho))
corr(max_Q_length', residuals_alpha_rho', 'Type', 'Spearman')

corrcoef(log(max_Q_length), log(residuals_beta))
corr(max_Q_length', residuals_beta', 'Type', 'Spearman')

figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');
ax = nexttile;
plot(log(residuals_alpha_patch), log(max_Q_length), 'o', 'markersize', 8, 'linewidth', 1)
title('alpha patch fit residuals vs max Q stay length (Extreme states)')
xlabel('alpha patch residuals')
ylabel('Length of extreme state')
set(ax,'xscale' ,'log')
set(ax,'yscale' ,'log')

ax = nexttile;
plot(log(residuals_alpha_rho), log(max_Q_length), 'o', 'markersize', 8, 'linewidth', 1)
title('alpha rho fit residuals vs max Q stay length (Extreme states)')
xlabel('alpha rho residuals')
ylabel('Length of extreme state')
set(ax,'xscale' ,'log')
set(ax,'yscale' ,'log')

ax = nexttile;
plot(log(residuals_beta), log(max_Q_length), 'o', 'markersize', 8, 'linewidth', 1)
title('beta fit residuals vs max Q stay length (Extreme states)')
xlabel('beta residuals')
ylabel('Length of extreme state')
set(ax,'xscale' ,'log')
set(ax,'yscale' ,'log')

%% checks - simulation and fitting scripts end up with the same outputs 
Block = 1;

TestParams = [0.4, 0.1, 5];

% simulate data
out = simulateM25_RLOn(TestParams, Block, Env); % change this depending on model to test
SimData.Action = out.Action;
SimData.PatchOrder = Env.PatchOrder{Block};
SimData.NextPatch = out.NextPatch; 
[SubjNLLEval, RecoveredData] = NLL_M25_RLOn(TestParams, Env, SimData);

% issue is next patch - Q values start to diverge as soon as there's a
% difference in the NextPatch prediction. I've now changed the next patch
% prediction to be just mode (not stochastic). If this improves recovery,
% then will need to log simulated next patch 