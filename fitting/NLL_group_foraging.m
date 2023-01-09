function [NLL, numObservations] = NLL_group_foraging(pars,Data, Patch_Order,environ, model_table, model)


numSubjects = numel(Data);      % number of subjects is the number of cells in the array
NLL = zeros(numSubjects,1);
numObservations = zeros(numSubjects,1);

for iS = 1:numSubjects
    if strcmp(model_table.policy{model}, 'MVT_learning')
        [NLL(iS), ~, numObservations(iS)] = NLL_subject_foraging_MVT(pars,Data{iS},environ);
    elseif strcmp(model_table.policy{model}, 'on')
        [NLL(iS), ~, ~, ~, ~, numObservations(iS)] = NLL_subject_foraging_RL_on(pars,Data{iS},Patch_Order{iS},environ,model_table, model);
    elseif strcmp(model_table.policy{model}, 'off')
        [NLL(iS), ~, ~, ~, ~, numObservations(iS)] = NLL_subject_foraging_RL_off(pars,Data{iS},Patch_Order{iS},environ, model_table, model);
    end
end