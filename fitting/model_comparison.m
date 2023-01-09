% model comparison of foraging
% Emma Scholeyl, 28 Nov 2022
clearvars
close all;
addpath(genpath('~/Dropbox/foraging-project/code'))
model_table = readtable('~/Dropbox/foraging-project/details/input_model_table.xlsx');

model = [2 21 25];
block = 2; % rich = 1, poor = 2
block_names = {'rich', 'poor'};
numParams  = 3;

% fmincon options
lowerBounds = [0,0,0];  % [alpha_Q, alpha_rho, beta] % parameter bounds
upperBounds = [1,1,100];  % arbitrary upper bound on beta to stop pathological behaviour
%options = optimoptions('fmincon','Display','none', 'Algorithm', 'active-set'); % don't display
options = optimoptions('fmincon','Display','none'); % don't display

%% model recovery - plot confusion matrix
% extent to which simulated data from model A can be explained by model B

% set up the foraging environment
set_up_environ % get the environment parameters
if block == 1
    patch_order = environ.patch_order.rich;
elseif block == 2
    patch_order = environ.patch_order.poor;
end

CM = zeros(size(model,2));

for count = 1:100
    count

    figure(1); clf;
    FM = round(100*CM/sum(CM(1,:)))/100;
    t = imageTextMatrix(FM);
    set(t(FM'<0.3), 'color', 'w')
    hold on;
    [l1, l2] = addFacetLines(CM);
    set(t, 'fontsize', 22)
    title(['count = ' num2str(count)]);
    set(gca, 'xtick', [1:5], 'ytick', [1:5], 'fontsize', 28, ...
        'xaxislocation', 'top', 'tickdir', 'out')
    xlabel('fit model')
    ylabel('simulated model')

    drawnow

    for m = 1:size(model, 2) % for each model we want to test
        sim_model = model(m);
        testparams = [rand, 0.0001+rand(1,1)*(0.2-0.0001), 1+exprnd(2)]; % randomly select parameters each time

        % simulate data
        if strcmp(model_table.policy{sim_model}, 'MVT_learning')
            [output, nObv] = simulate_MVT(testparams,environ,block,model_table,sim_model);

        elseif strcmp(model_table.policy{sim_model}, 'on')
            [output, ~, ~, ~, nObv]  = simulate_RL_on_policy(testparams,environ,block,model_table,sim_model);

        elseif strcmp(model_table.policy{sim_model}, 'off')
            [output, ~, ~, ~, nObv]  = simulate_RL_off_policy(testparams,environ,block,model_table,sim_model);

        end
        sim_foraging = [output.time_in_block, output.a_selected, output.patch_reward];

        for mi = 1:size(model,2)
            fit_model = model(mi);
            % fit all of the models back to the simulated data of each model
            if strcmp(model_table.policy{fit_model}, 'MVT_learning')
                f = @(x)NLL_subject_foraging_MVT(x,sim_foraging,environ);
            elseif strcmp(model_table.policy{fit_model}, 'on')
                f = @(x)NLL_subject_foraging_RL_on(x,sim_foraging,patch_order,environ,model_table,fit_model);
            elseif strcmp(model_table.policy{fit_model}, 'off')
                f = @(x)NLL_subject_foraging_RL_off(x,sim_foraging,patch_order,environ,model_table,fit_model);
            end

            params0 = [rand, 0.0001+rand(1,1)*(0.2-0.0001), 1+exprnd(2)];
            [~,NLL_eval] = fmincon(f,params0,[],[],[],[],lowerBounds,upperBounds, [], options); % start search from random location each time
            BIC(mi) = numParams * log(nObv) + 2*NLL_eval;
            [M, iBEST] = min(BIC);
            BEST = BIC == M;
            BEST = BEST / sum(BEST);
        end
                    clear BIC M iBEST

        CM(m,:) = CM(m,:) + BEST;
    end
end

%
figure(1); 
title('')
set(gcf, 'Position', [811   417   500   400])
set(gca, 'fontsize', 28);


%% AIC/BIC
models_AIC = zeros(size(model));
models_BIC = zeros(size(model));

for m = 1:size(model, 2)
    load(sprintf('~/Dropbox/foraging-project/results/M%d/fitting_results_%s', model(m), block_names{block}));
    sum_AIC = sum(AIC);
    sum_BIC = sum(BIC);
    models_AIC(m) = sum_AIC;
    models_BIC(m) = sum_BIC;
end

% plot

figure
bar(models_BIC)
xticklabels(model_table.policy(model))
ylabel('BIC')

figure
bar(models_AIC)
xticklabels(model_table.policy(model))
ylabel('AIC')


