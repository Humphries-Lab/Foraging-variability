function  [stats, sum_reward, mean_Q_stay, mean_Q_leave] = extract_stats(data, model, block)

% FOR EACH SIMULATION

% index to extract rows where agent is leaving
leaving_index = logical([diff(data.a_selected)==-1; 0]); % find the stay action where it changes to leave in the next state

% extract leaving times and the patch type
lt = [data.time_in_patch(leaving_index), data.patch_type(leaving_index)];
[lt_mean, lt_sd, lt_grps] = grpstats(lt(:,1), lt(:,2), ["mean", "std", "gname"]);
lt_grps = str2num(cell2mat(lt_grps));

% extract leaving reward rate threshold
lrr = [data.patchRR(leaving_index), data.patch_type(leaving_index)];
[lrr_mean, lrr_sd] = grpstats(lrr(:,1), lrr(:,2), ["mean", "std"]);

% extract average reward rate estimate at leaving (leaving rho, lrho)
lrho = [data.rho(leaving_index), data.patch_type(leaving_index)];
[lrho_mean, lrho_sd] = grpstats(lrho(:,1), lrho(:,2), ["mean", "std"]);

% correct if bugs out
if size(lt_mean,1) ~= 3 
    lt_mean = [nan; nan; nan]; lt_sd = [nan; nan; nan]; lrr_mean = [nan; nan; nan]; lrr_sd = [nan; nan; nan]; lrho_mean = [nan; nan; nan]; lrho_sd = [nan; nan; nan]; lt_grps = [nan; nan; nan];
end

stats = [lt_mean, lt_sd, lrr_mean, lrr_sd, lrho_mean, lrho_sd, lt_grps, block*ones(3,1), model*ones(3,1)]; % add column defining model type

stats = array2table(stats, 'VariableNames', {'lt_mean', 'lt_sd', 'lrr_mean', 'lrr_sd', 'lrho_mean', 'lrho_sd', 'patch_type', 'environ_type', 'model'});

% how much reward did they get?
sum_reward = sum(data.patch_reward);

%extract Q for each patch if R-learning model
if ~ismember(model, [1, 2]) % if not an MVT model
    for i = 1:max(data.patch_number) % for each patch
        Q_stay_i = (data.a_selected == 2 & data.patch_number == i); % index
        Q.stay{1+i} = [data.Q_stay(Q_stay_i), data.patch_type(Q_stay_i), data.time_in_patch(Q_stay_i)+1];

        Q_leave_i = (data.a_selected == 1 & data.patch_number == i); % index
        Q.leave{1+i} = [data.Q_leave(Q_leave_i), data.patch_type(Q_leave_i), (1:size(data.Q_leave(Q_leave_i)))'];
    end
    Q_stay_all = cat(1, Q.stay{:});
    [mean_Q_stay, grps] = grpstats(Q_stay_all(:,1), [Q_stay_all(:,2), Q_stay_all(:,3)], ["mean", "gname"]);
    mean_Q_stay = [mean_Q_stay, str2double(grps), block*ones(length(mean_Q_stay), 1), model*ones(length(mean_Q_stay), 1)];
    mean_Q_stay = array2table(mean_Q_stay, 'VariableNames', {'Q_stay', 'patch_type', 'state', 'environ_type', 'model'});

    Q_leave_all = cat(1, Q.leave{:});
    [mean_Q_leave, grps] = grpstats(Q_leave_all(:,1), [Q_leave_all(:,2), Q_leave_all(:,3)], ["mean", "gname"]);
    mean_Q_leave = [mean_Q_leave, str2double(grps), block*ones(length(mean_Q_leave), 1), model*ones(length(mean_Q_leave), 1)];
    mean_Q_leave = array2table(mean_Q_leave, 'VariableNames', {'Q_leave', 'patch_type', 'state', 'environ_type', 'model'});

end
end

