% fitting RL model functions to foraging data
% Emma Scholey 9 Jun 2022
clearvars
close all;
addpath(genpath('~/Dropbox/foraging-project/code'))
model_table = readtable('~/Dropbox/foraging-project/details/input_model_table.xlsx');

model = 25;

set_up_environ % get the environment parameters
block = 1; % rich = 1, poor = 2

% fmiamescon options
lowerBounds = [0,0,0];  % [alpha_Q, alpha_rho, beta] % parameter bounds
upperBounds = [1,1,100];  % arbitrary upper bound on beta to stop pathological behaviour
%options = optimoptions('fmincon','Display','none', 'Algorithm', 'active-set'); % don't display
options = optimoptions('fmincon','Display','none'); % don't display

% get their data
subj = [1:23 25:40];

[Data, Patch_Order] = prepare_subject_data(subj, block, environ, r);

numSubjects = size(Data,2);
numParams = length(lowerBounds);

 %% test NLL computation for subject
% % fit one subject for testing
% subj = 5;
% params0 = [0.7 0.3 10];
% 
% % test NLL is working
% if strcmp(model_table.policy{model}, 'MVT_learning')
%     [NLL_eval] = NLL_subject_foraging_MVT(params0,Data{subj},environ);
% elseif strcmp(model_table.policy{model}, 'on')
%     [NegLogLikelihood,rho,Q_stay,Q_leave,PAction,numObservations] = NLL_subject_foraging_RL_on(params0,Data{subj},Patch_Order{subj},environ, model_table, model);
% elseif strcmp(model_table.policy{model}, 'off')
%     [NegLogLikelihood,rho,Q_stay,Q_leave,PAction,numObservations] = NLL_subject_foraging_RL_off(params0,Data{subj},Patch_Order{subj},environ, model_table, model);
% end

%% grid search for subjects - marginal likelilood distributions

subj = [1:23 25:40];

% range of values for grid search of starting points for search
alpha_vector = 0:0.05:1;
beta_vector = 0:2.5:50;

% create grid on which to initialise searches
[alpha_patch_grid,alpha_rho_grid, beta_grid] = ndgrid(alpha_vector,alpha_vector,beta_vector);

% create storage
N_alpha = numel(alpha_vector);
N_beta = numel(beta_vector);

NLL_eval = zeros(numSubjects,N_alpha,N_alpha,N_beta);

for iAlphaPatch = 1:N_alpha
    for iAlphaRho = 1:N_alpha
        for iBeta = 1:N_beta
            params0 = [alpha_patch_grid(iAlphaPatch,iAlphaRho,iBeta),alpha_rho_grid(iAlphaPatch,iAlphaRho,iBeta), beta_grid(iAlphaPatch,iAlphaRho,iBeta)]  % get the current set of initial search parameters
            for iData = 1:numSubjects
                if strcmp(model_table.policy{model}, 'MVT_learning')
                    [NLL_eval(iData, iAlphaPatch,iAlphaRho,iBeta)] = NLL_subject_foraging_MVT(params0,Data{iData},environ);
                elseif strcmp(model_table.policy{model}, 'on')
                    [NLL_eval(iData, iAlphaPatch,iAlphaRho,iBeta)] = NLL_subject_foraging_RL_on(params0,Data{iData},Patch_Order{iData},environ, model_table, model);
                elseif strcmp(model_table.policy{model}, 'off')
                    [NLL_eval(iData, iAlphaPatch,iAlphaRho,iBeta)] = NLL_subject_foraging_RL_off(params0,Data{iData},Patch_Order{iData},environ, model_table, model);
                end
            end
        end
    end
end

%% find min NLL and plot ML
for iData = 1:numSubjects
    % for each subject, find their minimum NLL across grid of starting points
    their_NLL = squeeze(NLL_eval(iData,:,:,:));
    minNLL(iData) = min(min(min(their_NLL)));   % minimum negative log likelihood over all starting positions
    [ix_alpha_patch,ix_alpha_rho,ix_beta] = ind2sub(size(alpha_patch_grid),find(their_NLL == minNLL(iData)));    % indices of location of minimum
    ix_alpha_patch(iData) = ix_alpha_patch(1); %Note: if more than one starting location converges on same parameters then will have same NLL, so use first
    ix_alpha_rho(iData) = ix_alpha_rho(1);
    ix_beta(iData) = ix_beta(1);

    % get the corresponding fitted parameter values.
    minNLL_params(iData,:) = [alpha_patch_grid(ix_alpha_patch(iData),ix_alpha_rho(iData),ix_beta(iData)), alpha_rho_grid(ix_alpha_patch(iData), ix_alpha_rho(iData),ix_beta(iData)),beta_grid(ix_alpha_patch(iData), ix_alpha_rho(iData),ix_beta(iData))];
end

% plot parameter distributions - marginal likelihoods
lik = exp(-NLL_eval); % compute the likelihood, (rather than the log (negative LL))

close all
figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');
for iData = 1:4
    nexttile;
    their_lik = squeeze(lik(iData,:,:,:));
    patch_tmp = sum(sum(their_lik, 3),2); % sum ml for alpha_patch across all beta and alpha_rho (MARGINAL likelihood)
    alpha_patch_marglik = patch_tmp/sum(patch_tmp);
    rho_tmp = sum(sum(their_lik, 3),1); % do the same for alpha_rho
    alpha_rho_marglik = rho_tmp/sum(rho_tmp);
    plot(alpha_vector, alpha_patch_marglik,alpha_vector, alpha_rho_marglik)
    xlabel('alpha', 'FontSize', 12)
    ylabel('p(data|model)', 'FontSize', 12)
    ylim([0 1])
    title(sprintf('Subject %d', iData))

end
leg = legend({'Patch learning rate', 'Rho learning rate'}, 'FontSize', 14);
leg.Layout.Tile = 'north';
title(tl, 'marginal likelihood distribution - learning rates', 'FontSize', 16)

figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');

for iData = 1:4
    nexttile;
    their_lik = squeeze(lik(iData,:,:,:));
    beta_tmp = sum(squeeze(sum(their_lik, 1)),1); % do the same for betas
    beta_marglik = beta_tmp/sum(beta_tmp);
    plot(beta_vector, beta_marglik)
    xlabel('beta', 'FontSize', 12)
    ylim([0 1])
    ylabel('p(data|model)', 'FontSize', 12)
    title(sprintf('Subject %d', iData))

end
leg = legend({'Softmax temperature' }, 'FontSize', 14);
leg.Layout.Tile = 'north';
title(tl, 'marginal likelihood distribution - softmax temperature', 'FontSize', 16)

clear ix_beta ix_alpha_patch ix_alpha_rho


%% test NLL fitting (fmincon) for one subject

% params0 = [0.8500    0.6500    1.0000]; % [alpha_Q, alpha_rho, beta]
% % test NLL is working
% if strcmp(model_table.policy{model}, 'MVT_learning')
%     f = @(x)NLL_subject_foraging_MVT(x,Data{subj},environ);
% elseif strcmp(model_table.policy{model}, 'on')
%     f = @(x)NLL_subject_foraging_RL_on(x,Data{subj},Patch_Order{subj},environ, model_table, model);
% elseif strcmp(model_table.policy{model}, 'off')
%     f = @(x)NLL_subject_foraging_RL_off(x,Data{subj},Patch_Order{subj},environ,model_table, model);
% end
% 
% [params_fitted, NLL_eval, ~, ~, ~, ~, hessian] = fmincon(f,params0,[],[],[],[],lowerBounds,upperBounds, [], options);
% 
% SE_params_fitted = sqrt(inv(hessian)); % The diagonal terms of H−1 correspond to variances for each parameter separately, and their square roots measure one standard error on the parameter. - Daw 2011

%% test few starting points for one subject
% alpha_vector = 0.05:0.2:1;
% beta_vector = 1:2:11;
%
% % create grid on which to initialise searches
% [alpha_patch_grid,alpha_rho_grid, beta_grid] = ndgrid(alpha_vector,alpha_vector,beta_vector);
%
% % create storage
% N_alpha = numel(alpha_vector);
% N_beta = numel(beta_vector);
%
% NLL_eval = zeros(N_alpha,N_alpha,N_beta);
% params_fitted = zeros(N_alpha,numParams,N_alpha,N_beta);
%
% for iAlphaPatch = 1:N_alpha
%     for iAlphaRho = 1:N_alpha
%         for iBeta = 1:N_beta
%             params0 = [alpha_patch_grid(iAlphaPatch,iAlphaRho,iBeta),alpha_rho_grid(iAlphaPatch,iAlphaRho,iBeta), beta_grid(iAlphaPatch,iAlphaRho,iBeta)]  % get the current set of initial search parameters
%             if strcmp(model_table.policy{model}, 'MVT_learning')
%                 f = @(x)NLL_subject_foraging_MVT(x,Data{subj},environ);
%             elseif strcmp(model_table.policy{model}, 'on')
%                 f = @(x)NLL_subject_foraging_RL_on(x,Data{subj},Patch_Order{subj},environ,model_table, model);
%             elseif strcmp(model_table.policy{model}, 'off')
%                 f = @(x)NLL_subject_foraging_RL_off(x,Data{subj},Patch_Order{subj},environ,model_table, model);
%             end
%             [params_fitted(iAlphaPatch, :, iAlphaRho, iBeta), NLL_eval(iAlphaPatch, iAlphaRho, iBeta)] = fmincon(f,params0,[],[],[],[],lowerBounds,upperBounds, [], options);
%         end
%     end
% end
%
% minNLL = min(NLL_eval(:));
% % get indices of their min NLL
% [ix_alpha_patch, ix_alpha_rho,ix_beta] = ind2sub(size(alpha_patch_grid),find(NLL_eval == minNLL));    % indices of location of minimum. Do it this way because just using find doesn't work with 3D matrix
% % get the corresponding fitted parameter values.
% minNLL_params = squeeze(params_fitted(ix_alpha_patch,:,ix_alpha_rho,ix_beta));
%
% % parameter distribution for single subject
% close all
% figure
% histogram(params_fitted(:,1,:,:)); % alpha patch
% xlabel('alpha patch - fitted')
% title('distribution of fit parameter values - alpha patch')
%
% figure
% histogram(params_fitted(:,2,:,:));  % alpha rho
% xlabel('alpha rho - fitted')
% title('distribution of fit parameter values - alpha rho')
%
% figure
% histogram(params_fitted(:,3,:,:)); % beta
% xlabel('beta - fitted')
% title('distribution of fit parameter values - beta')

%% testing group parameter fitting
% subj = [1:23 25:40];
%
% params0 = [0.26 0.05 8]; % [alpha_Q, alpha_rho, beta]
%
% % test NLL computation for all subjects
% [NLL] = NLL_group_foraging(params,Data, Patch_Order,environ, model_table, model);
% group_NLL = sum(NLL);

% test fitting to whole group
% f = @(x)NLL_group_foraging(x,Data, Patch_Order,environ, model_table, model);
% [group_params_fitted,group_NLL_eval] = fmincon(f,params0,[],[],[],[],lowerBounds,upperBounds, options);

%% test fitting to each subject of group
% params0 = [0.3, 0.1, 4];  % initial search parameters
% each_params_fitted = zeros(numSubjects,numParams);
% each_NLL_eval = zeros(numSubjects,1);
%
% for iS = 1:numSubjects
%     if strcmp(model_table.policy{model}, 'MVT_learning')
%         f = @(x)NLL_subject_foraging_MVT(x,Data{iS},environ);
%     elseif strcmp(model_table.policy{model}, 'on')
%         f = @(x)NLL_subject_foraging_RL_on(x,Data{iS},Patch_Order{iS}, environ,model_table, model);
%     elseif strcmp(model_table.policy{model}, 'off')
%         f = @(x)NLL_subject_foraging_RL_off(x,Data{iS},Patch_Order{iS}, environ,model_table, model);
%     end
%     [each_params_fitted(iS,:),each_NLL_eval(iS)] = fmincon(f,params0,[],[],[],[],lowerBounds,upperBounds, options);
% end

%% fitting for each person in group with different starting points
nStarts = 50; % how many starting points to avoid local minima

each_params_fitted = zeros(numSubjects,numParams,nStarts);
each_NLL_eval = zeros(numSubjects,nStarts);

for ii = 1:nStarts
    params0(ii,:) = [rand, rand, exprnd(2)] % randomly select parameters each time
    for iData = 1:numSubjects
        if strcmp(model_table.policy{model}, 'MVT_learning')
            f = @(x)NLL_subject_foraging_MVT(x,Data{iData}, environ);
        elseif strcmp(model_table.policy{model}, 'on')
            f = @(x)NLL_subject_foraging_RL_on(x,Data{iData},Patch_Order{iData}, environ,model_table, model);
        elseif strcmp(model_table.policy{model}, 'off')
            f = @(x)NLL_subject_foraging_RL_off(x,Data{iData},Patch_Order{iData}, environ,model_table, model);
        end
        [each_params_fitted(iData,:,ii),each_NLL_eval(iData,ii), ~, ~, ~, ~, hessian(:, :, iData, ii)] = ...
            fmincon(f,params0(ii,:),[],[],[],[],lowerBounds,upperBounds, [], options);

        % calculate standard error on the parameter fits using the hessian
        % matrix
        tmp = sqrt(diag(inv(hessian(:,:,iData,ii))))'; % The diagonal terms of H−1 correspond to variances for each parameter separately, and their square roots measure one standard error on the parameter. - Daw 2011
        tmp(~isreal(tmp),:) = nan;
        SE_params_fitted(iData,:,ii) = tmp;
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
[~, numObservations] = NLL_group_foraging([0 0 1],Data, Patch_Order, environ, model_table, model);

for iData = 1:numSubjects
    % for each subject, find their minimum NLL across grid of starting points
    their_NLL = squeeze(each_NLL_eval(iData,:));
    minNLL(iData) = min(their_NLL);   % minimum negative log likelihood over all starting positions
    ix = find(their_NLL == minNLL(iData));    % indices of location of minimum

    % get the corresponding fitted parameter values. Note: if more than one starting location converges on same parameters then will have same NLL, so use first
    minNLL_params_fitted(iData,:) = squeeze(each_params_fitted(iData,:,ix));

    % get the corresponding SE for the fitted parameter.
    minNLL_SE_params_fitted(iData,:) = squeeze(SE_params_fitted(iData,:,ix));

    % BIC for each person
    BIC(iData) = numParams * log(numObservations(iData)) + 2*minNLL(iData);
    AIC(iData) = 2/numObservations(iData) * minNLL(iData) + 2 * numParams/numObservations(iData);
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
for iData = 1:4
    ax = nexttile;
    histogram(squeeze(each_params_fitted(iData,1,:))); % alpha patch
    title(sprintf('Subject %d', iData))
    xlabel('alpha patch - fitted', 'FontSize', 12)
end
title(tl, 'distribution of fit parameter values - alpha patch', 'FontSize', 16);

figure;
tl = tiledlayout('flow', 'TileSpacing', 'Compact');
for iData = 1:4
    nexttile;
    histogram(each_params_fitted(iData,2,:), 'FaceColor',[0.8500 0.3250 0.0980]);  % alpha rho
    title(sprintf('Subject %d', iData))
    xlabel('alpha rho - fitted', 'FontSize', 12)
end
title(tl, 'distribution of fit parameter values - alpha rho', 'FontSize', 16);

figure;
tl = tiledlayout('flow', 'TileSpacing', 'Compact');
for iData = 1:4
    nexttile;
    histogram(each_params_fitted(iData,3,:)); % beta
    title(sprintf('Subject %d', iData))
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

for iData = 1:4 % just test on a few participants
    nexttile;
    scatter3(squeeze(each_params_fitted(iData,1,:)),squeeze(each_params_fitted(iData,2,:)),squeeze(each_params_fitted(iData,3,:)), 50, squeeze(each_params_fitted(iData,3,:)), 'filled');
    hold on;
    scatter3(params0(:,1),params0(:,2),params0(:,3), 50, params0(:,3), "r", 'filled');
    colormap("winter");
    colorbar;
    title(sprintf('Subject %d', iData))
    xlabel('alpha patch', 'FontSize', 12)
    ylabel('alpha rho', 'FontSize', 12)
    zlabel('beta', 'FontSize', 12)
    % plot lines between the start and fit points
    for ii = 1:nStarts
        v1 = [params0(ii,1),params0(ii,2),params0(ii,3)];
        v2 = [squeeze(each_params_fitted(iData,1,ii)),squeeze(each_params_fitted(iData,2,ii)),squeeze(each_params_fitted(iData,3,ii))];
        v = [v2;v1];
        plot3(v(:,1), v(:,2),v(:,3),'r', 'LineStyle', '--')
    end
end
title(tl, 'Distance from start to fit parameters', 'FontSize', 16)
leg = legend({'Fit parameters', 'Start parameters'}, 'FontSize', 14);
leg.Layout.Tile = 'north';

%% parameter recovery
% match the range of your simulation parameters to the range of values obtained by your fit
% if block == 1
%     patch_order = environ.patch_order.rich;
% elseif block == 2
%     patch_order = environ.patch_order.poor;
% end
%
% testparams = [0.4 0.005 8]; % alpha_patch, alpha_rho, beta
%
% if strcmp(model_table.policy{model}, 'MVT_learning')
%     output = simulate_MVT(testparams,environ,block,model_table,model);
% elseif strcmp(model_table.policy{model}, 'on')
%     output = simulate_RL_on_policy(testparams,environ,block,model_table,model);
% elseif strcmp(model_table.policy{model}, 'off')
%     output = simulate_RL_off_policy(testparams,environ,block,model_table,model);
% end
%
% sim_foraging = [output.time_in_block, output.a_selected, output.patch_reward];
%
% % recover data
% for ii = 1:100 % for 100 random starting parameters
%     params0 = [rand, rand, exprnd(3)];
%     % recover data
%     if strcmp(model_table.policy{model}, 'MVT_learning')
%         f = @(x)NLL_subject_foraging_MVT(x,sim_foraging,environ);
%     elseif strcmp(model_table.policy{model}, 'on')
%         f = @(x)NLL_subject_foraging_RL_on(x,sim_foraging,patch_order,environ,model_table,model);
%     elseif strcmp(model_table.policy{model}, 'off')
%         f = @(x)NLL_subject_foraging_RL_off(x,sim_foraging,patch_order,environ,model_table,model);
%     end
%
%     [sim_params_fitted,NLL_eval] = fmincon(f,params0,[],[],[],[],lowerBounds,upperBounds, [], options); % start search from random location each time
%
%     fit_params(1,ii) = sim_params_fitted(1); % alpha patch rr fit
%     fit_params(2,ii) = sim_params_fitted(2); % alpha rho fit
%     fit_params(3,ii) = sim_params_fitted(3); % beta fit
% end

%% parameter recovery for range of values
if block == 1
    patch_order = environ.patch_order.rich;
elseif block == 2
    patch_order = environ.patch_order.poor;
end

n_sim = 250; % number of iterations

sim_params = zeros([numParams,n_sim]);
fit_params = zeros([numParams,n_sim]);

for ii = 1:n_sim
    testparams = [0.001+rand(1,1)*(0.6-0.001), 0.001+rand(1,1)*(0.15-0.001), exprnd(10)] % randomly select parameters each time

    % simulate data
    if strcmp(model_table.policy{model}, 'MVT_learning')
        [output, nObv(ii)] = simulate_MVT(testparams,environ,block,model_table,model);
    elseif strcmp(model_table.policy{model}, 'on')
        [output, Q_stay{ii}, ~, ~, nObv(ii)] = simulate_RL_on_policy(testparams,environ,block,model_table,model);
    elseif strcmp(model_table.policy{model}, 'off')
        [output, Q_stay{ii}, ~, ~, nObv(ii)] = simulate_RL_off_policy(testparams,environ,block,model_table,model);
    end
    sim_foraging = [output.time_in_block, output.a_selected, output.patch_reward];

    % recover data
    if strcmp(model_table.policy{model}, 'MVT_learning')
        f = @(x)NLL_subject_foraging_MVT(x,sim_foraging,environ);
    elseif strcmp(model_table.policy{model}, 'on')
        f = @(x)NLL_subject_foraging_RL_on(x,sim_foraging,patch_order,environ,model_table,model);
    elseif strcmp(model_table.policy{model}, 'off')
        f = @(x)NLL_subject_foraging_RL_off(x,sim_foraging,patch_order,environ,model_table,model);
    end

    for nStarts = 1:10 % do a few starting points to try and avoid global minima
        params0 = [0.001+rand(1,1)*(0.6-0.001), 0.001+rand(1,1)*(0.15-0.001), exprnd(10)];
        [sim_params_fitted(nStarts,:),NLL_eval(nStarts)] = fmincon(f,params0,[],[],[],[],lowerBounds,upperBounds, [], options); % start search from random location each time
    end

    % get the global minima
    ix = find(NLL_eval == min(NLL_eval));    % indices of location of minimum negative log likelihood over all starting positions
    % get the corresponding fitted parameter values. Note: if more than one starting location converges on same parameters then will have same NLL, so use first
    minNLL_sim_params_fitted = sim_params_fitted(ix,:);

    sim_params(1,ii) = testparams(1); % alpha patch rr simulated
    sim_params(2,ii) = testparams(2); % alpha rho simualted
    sim_params(3,ii) = testparams(3); % beta simulated

    fit_params(1,ii) = minNLL_sim_params_fitted(1); % alpha patch rr fit
    fit_params(2,ii) = minNLL_sim_params_fitted(2); % alpha rho fit
    fit_params(3,ii) = minNLL_sim_params_fitted(3); % beta fit

    % get their Q_stay table (based on the fit parameter values)
    if strcmp(model_table.policy{model}, 'on')
        [~,~,Q_stay_tables{ii},~,~,~] = NLL_subject_foraging_RL_on(minNLL_sim_params_fitted,sim_foraging,patch_order,environ,model_table,model);
    elseif strcmp(model_table.policy{model}, 'off')
        [~,~,Q_stay_tables{ii},~,~,~] = NLL_subject_foraging_RL_off(minNLL_sim_params_fitted,sim_foraging,patch_order,environ,model_table,model);
    end

    clear minNLL_sim_params_fitted sim_params_fitted NLL_eval
end

%% plot correlation - simulated value versus actual value
close all
names = {'learning rate - Q/patch' 'learning rate - rho' 'softmax temperature'};
symbols = {'\alpha' '\alpha' '\beta'};
figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');

% for each parameter, plot sim vs fit 
for i= 1:size(sim_params,1)
    ax = nexttile;
    plot(sim_params(i,:), fit_params(i,:), 'o', 'markersize', 8, 'linewidth', 1)
    corrcoef(sim_params(i,:),fit_params(i,:))
    corr(sim_params(i,:)',fit_params(i,:)', 'Type', 'Spearman')
    title(sprintf('%s', names{i}))
    xlabel(sprintf('simulated %s', symbols{i}))
    ylabel(sprintf('fit %s', symbols{i}))
end
title(tl, 'Parameter recovery - correlations between simulated and fit parameters')
set(ax,'xscale', 'log', 'yscale' ,'log')

corrcoef(log(sim_params(3,:)),log(fit_params(3,:)))

figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');
nexttile;
plot(fit_params(1,:), fit_params(2,:), 'o', 'markersize', 8, 'linewidth', 1)
title('fit alpha Q/patch vs fit alpha rho')
xlabel('fit alpha Q/patch')
ylabel('fit alpha rho')

ax = nexttile;
plot(fit_params(2,:), fit_params(3,:), 'o', 'markersize', 8, 'linewidth', 1)
title('fit alpha rho vs fit beta')
xlabel('fit alpha rho')
ylabel('fit beta')
set(ax,'yscale' ,'log')

ax = nexttile;
plot(fit_params(1,:), fit_params(3,:), 'o', 'markersize', 8, 'linewidth', 1)
title('fit alpha Q/patch vs fit beta')
xlabel('fit alpha Q/patch')
ylabel('fit beta')
set(ax,'yscale' ,'log')
title(tl, 'Parameter trade-offs - correlations between fit parameters')

%% Checks for parameter recovery - extreme states 
for ii = 1:n_sim
    max_Q_length(ii) = max(sum(Q_stay{ii} ~= 0));
end

residuals_alpha_patch = sqrt((fit_params(1,:) - sim_params(1,:)).^2);
residuals_alpha_rho = sqrt((fit_params(2,:) - sim_params(2,:)).^2);
residuals_beta = sqrt((fit_params(3,:) - sim_params(3,:)).^2);

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
if block == 1
    patch_order = environ.patch_order.rich;
elseif block == 2
    patch_order = environ.patch_order.poor;
end

testparams = [0.7, 0.3, 20];

% simulate data
if strcmp(model_table.policy{model}, 'on')
    [output, Q_stay, Q_leave, ~, nObv] = simulate_RL_on_policy(testparams,environ,block,model_table,model);
elseif strcmp(model_table.policy{model}, 'off')
    [output, Q_stay, ~, ~, nObv] = simulate_RL_off_policy(testparams,environ,block,model_table,model);
end
sim_foraging = [output.time_in_block, output.a_selected, output.patch_reward];

% get their Q_stay table (based on the fit parameter values)
if strcmp(model_table.policy{model}, 'on')
    [~,~,Q_stay_tables,Q_leave_fit,~,nObv_fit] = NLL_subject_foraging_RL_on(testparams,sim_foraging,patch_order,environ,model_table,model);
elseif strcmp(model_table.policy{model}, 'off')
    [~,~,Q_stay_tables,Q_leave_fit,~,nObv_fit] = NLL_subject_foraging_RL_off(testparams,sim_foraging,patch_order,environ,model_table,model);
end

% only issue is one extra step in the simulated data. Can correct that
% easily, and don't think that's causing the issues with recovery. 