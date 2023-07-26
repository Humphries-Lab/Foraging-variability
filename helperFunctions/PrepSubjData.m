function Data = PrepSubjData(SubjLT, block, Env, blockFlag)

%% build the foraging table from their MAIN session
% specify actions
Leave = 1;
Stay = 2;

if contains(blockFlag, 'combined') % if fitting as a continuous task (no separate blocks) 
    LT = SubjLT.leaveT; % pull out leaving times - note that this will do it in the correct block order for the participant [1 2] or [2 1]
    PatchOrder = SubjLT.patch; % pull out patch order
elseif contains(blockFlag, 'separate')
    LT = SubjLT.leaveT(SubjLT.env == block); % pull out leaving times
    PatchOrder = SubjLT.patch(SubjLT.env == block); % pull out patch order
end

LT = round(LT/Env.TimeStep)*Env.TimeStep; % round to nearest timestep (state precision)

% transform leaving times into stay/leave actions for each state
for ii = 1:length(LT)
    a{ii} = repelem([Stay Leave], round([LT(ii) Env.TravelTime]./Env.TimeStep));
end

A = cat(2, a{:})'; % concatenate all patches into one long block

% Block type
if contains(blockFlag, 'combined')
    B = SubjLT.env;
elseif contains(blockFlag, 'separate')
    B = SubjLT.env(SubjLT.env == block);
end

% calculate their experienced avgRR for each block (based on type & length in patch) 
% for iB = 1:2
%     blockLT = LT(SubjLT.env==iB); % only get LTs for correct block
%     blockPatchOrder = PatchOrder(SubjLT.env == iB);
%     for ii = 1:length(blockLT) % for each patch
%         reward_vector = Env.R(1:blockLT(ii),blockPatchOrder(ii));
%         patch_reward(ii) = sum(reward_vector); % sum reward for that patch, according to time in patch
%     end
% 
%     Data.experiencedAvgRR{iB} = sum(patch_reward)/(sum(blockLT)+Env.TravelTime*numel(blockLT)); % total reward/time
%     clear patch_reward
% end

Data.Action = A;
Data.PatchOrder = PatchOrder;
Data.NextPatch = 0; % this is required to force the NLL script to predict the next patch (for stochastic patch prediction models)
Data.BlockType = B; 
end