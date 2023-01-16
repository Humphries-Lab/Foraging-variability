
%% Script to run foraging simulations (MVT or R-learning model)
% Emma Scholey

%% Set up
clear
addpath(genpath('~/Dropbox/foraging/code'))

set_up_environ % generic script to set up foraging environment according to LeHeron et al (2020)

% model type - see table for model number to choose
save_data = 0;
n_sim = 100; % how many runs
num_models = 26;
model = 25;

model_table = readtable('~/Dropbox/foraging/details/input_model_table.xlsx');

T_behaviour = cell(num_models,1);
gT_behaviour = cell(num_models,1);
T_Q_stay = cell(num_models,1);
T_Q_leave = cell(num_models,1);
gT_Q_stay = cell(num_models,1);
gT_Q_leave = cell(num_models,1);

environ_type = {'rich', 'poor'};

%% run model
for model = model
    % tidy up first
    %close all
    vars = {'sim_stats','sim_data','sim_Q_stay','sim_Q_leave'};
    clear(vars{:})

    clear params
    % agent parameters - currently fitting rich and poor blocks separately
    params{1} = [0.4, 0.1, 5]; % [alpha_Q, alpha_rho, beta] - rich block
    params{2} = [0.4, 0.1, 5]; % [alpha_Q, alpha_rho, beta] - poor block

    %% optional - simulate real data for model validation

%     load(sprintf('~/Dropbox/foraging-project/results/M%d/fitting_results_rich.mat', model)) % load best fit parameters from real data
%     params{1} = minNLL_params_fitted;
%     load(sprintf('~/Dropbox/foraging-project/results/M%d/fitting_results_poor.mat', model)) % load best fit parameters from real data
%     params{2} = minNLL_params_fitted;
% 
%     n_sim = size(minNLL_params_fitted,1); % simulate each person

    %% Run
    for n = 1:n_sim
        n
        for block = 1:length(environ_type) % for each environment, [rich, poor]

            if contains(model_table.policy{model}, 'MVT')
                [sim_data.(environ_type{block}){n}] = simulate_MVT(params{block}(n,:), environ, block, model_table, model);
                [sim_stats.(environ_type{block}){n}, sum_reward.(environ_type{block}){n}] = extract_stats(sim_data.(environ_type{block}){n}, model, block);

            else % if R-learning models, also pull out Q-tables

                if strcmp(model_table.policy{model}, 'on')
                    % simulate main task
                    [sim_data.(environ_type{block}){n}, Q_stay_table.(environ_type{block}){n}, Q_leave_table.(environ_type{block}){n}] = ...
                        simulate_RL_on_policy(params{block}, environ, block, model_table, model);
                elseif strcmp(model_table.policy{model}, 'off')
                    % simulate main task
                    [sim_data.(environ_type{block}){n}, Q_stay_table.(environ_type{block}){n}, Q_leave_table.(environ_type{block}){n}] = ...
                        simulate_RL_off_policy(params{block}, environ, block, model_table, model);
                end

                [sim_stats.(environ_type{block}){n}, sum_reward.(environ_type{block}){n}, sim_Q_stay.(environ_type{block}){n}, sim_Q_leave.(environ_type{block}){n}] = ...
                    extract_stats(sim_data.(environ_type{block}){n}, model, block);
            end
        end
    end

    %% save runs and statistics to file
%     if save_data == 1
%         save(sprintf('~/Dropbox/foraging-project/data/sim_data/M%d.mat', model), 'sim_*')
%     end

    %% store leaving stats for all runs in table

    T_behaviour = [cat(1, sim_stats.rich{:}); cat(1, sim_stats.poor{:})];

    % group stats - across all runs
    gT_behaviour = grpstats(T_behaviour, ["environ_type", "patch_type"], ["mean", "meanci"], "DataVars", ["lt_mean", "lrr_mean", "lrho_mean"]);
    gT_behaviour.Properties.VariableNames = {...
        'environ_type',...
        'patch_type',...
        'GroupCount',...
        'g_mean_lt',...
        'g_ci_lt',...
        'g_mean_lrr',...
        'g_ci_lrr',...
        'g_mean_lrho',...
        'g_ci_lrho'};
    gT_behaviour.model = repmat(model, [6 1]); % add column defining model number

    %% Store Q stats
    if ~contains(model_table.policy{model}, 'MVT')
        T_Q_stay = [cat(1, sim_Q_stay.rich{:}); cat(1, sim_Q_stay.poor{:})];
        T_Q_leave = [cat(1, sim_Q_leave.rich{:}); cat(1, sim_Q_leave.poor{:})];

        gT_Q_stay = grpstats(T_Q_stay, ["environ_type", "patch_type", "state"], ["mean", "meanci"], "DataVars", "Q_stay");
        gT_Q_leave = grpstats(T_Q_leave, ["environ_type", "patch_type", "state"], ["mean", "meanci"], "DataVars", "Q_leave");

        gT_Q_stay.model = repmat(model, [size(gT_Q_stay, 1), 1]);
        gT_Q_leave.model = repmat(model, [size(gT_Q_leave, 1), 1]);

    end

    %% Store total reward stats
    T_reward = [cat(1,sum_reward.rich{:}), cat(1,sum_reward.poor{:})];
    mean_reward = mean(T_reward); % [rich poor]

    %% plot agent behaviour over single run
%     check_run = 30; % choose which run to plot
%     data = sim_data.poor{check_run};
%     %data = data(data.patch_number > 10, :);
%     plot_example_run(data, model)
%     %print(sprintf('~/Dropbox/foraging-project/results/M%d/M%d_example_run', model, model), '-dpng')

    %% store multiple model tests in one table
    all_models_mean_reward{model} = mean_reward;
    all_models_gT_behaviour{model} = gT_behaviour;
    all_models_T_behaviour{model} = T_behaviour;
    all_models_T_Q_stay{model} = T_Q_stay;
    all_models_T_Q_leave{model} = T_Q_leave;
    all_models_gT_Q_stay{model} = gT_Q_stay;
    all_models_gT_Q_leave{model} = gT_Q_leave;
end

%% save all model data together
% if save_data == 1
%     save('~/Dropbox/foraging-project/data/sim_data/all_models.mat','all_models_*')
% end
%% plot patch leaving times
load('~/Dropbox/foraging/raw_data/summary/young_variables/subMean_young.mat')
plot_leaving_stats_line(all_models_gT_behaviour{model}, subMean_young)