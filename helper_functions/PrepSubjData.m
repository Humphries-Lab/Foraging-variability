function Data = PrepSubjData(SubjLT, block, Env)

%% build the foraging table from their MAIN session
% specify actions
Leave = 1;
Stay = 2;

LT = SubjLT.leaveT(SubjLT.env == block); % pull out leaving times
LT = round(LT); % round to nearest integer (state precision)

PatchOrder = SubjLT.patch(SubjLT.env == block); % pull out patch order
tmp = repelem(PatchOrder, LT+Env.TravelTime); % create patch order for each state, for calculating reward later

% transform leaving times into stay/leave actions for each state
for ii = 1:length(LT)
    a{ii} = repelem([Stay Leave], [LT(ii) Env.TravelTime]);
end

A = cat(2, a{:})'; % concatenate all patches into one long block

Data.Action = A;
Data.PatchOrder = PatchOrder;
Data.NextPatch = 0; % this is required to force the NLL script to predict the next patch (for stochastic patch prediction models)
end