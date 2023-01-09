function [Data, patch_order] = prepare_subject_data(subj, env, environ, r)

addpath(genpath('~/Dropbox/foraging-project/data'))

%% build the foraging table from their MAIN session
load ~/Dropbox/foraging-project/data/raw_data/summary/young_variables/t_young.mat

group = {'Young_HC_ReDo/yc%g_forage.mat', subj};

for k = 1:length(group{2})
    % load derived participant data - leaving times for each subject by patch
    % and environ type

    subj_lt = t_young.leaveT(t_young.subj == k & t_young.env == env); % pull out leaving times for each subject
    patch_order{k} = t_young.patch(t_young.subj == k & t_young.env == env); % pull out patch type for each subject
    subj_lt(:,1) = round(subj_lt(:,1)); % round to nearest integer (state precision)

    T = {};

    % create time_in_patch column based on subject leaving times
    for ii = 1:length(subj_lt)
        lt = subj_lt(ii);
        T{ii} = [1:lt, repelem(NaN, environ.travel_time)];
    end
    T = cat(2, T{:})';

    % create patch type column. Basing this on their actual measured patch type, rather than the original master list, because each subject's patch order may not exactly
    % match the original (due to excluded patches where they 'fell asleep')
    patch_type = repelem(patch_order{k}, subj_lt+6); % convert patch type into states (e.g. what's happening at each timestep)

    % create actions column, based on subjects' leaving times
    A = ~isnan(T)+1; % convert patch/travel times to actions (stay = 2, leave = 1)

    % create reward column, based on their leaving times, patch type and time in patch
    R = zeros([length(T), 1]);
    for ii = 1:length(T)
        if isnan(T(ii))
            R(ii) = 0;
        else
            R(ii) = r(T(ii), patch_type(ii));
        end
    end

    % put everything together
    Data{k} = [T, A, R, patch_type];
end