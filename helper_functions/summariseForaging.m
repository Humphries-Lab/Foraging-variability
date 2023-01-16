function  [LTMean, LRRMean, TotalReward] = summariseForaging(run)

% FOR EACH SIMULATION

% extract leaving times for each patch type
for ii = 1:3
    LTMean(ii) = mean(run.LeavingTime(run.PatchOrder == ii));
    LTStandardDev(ii) = std(run.LeavingTime(run.PatchOrder == ii));
    LRRMean(ii) = mean(run.LeavingRR(run.PatchOrder == ii));
    LRRStandardDev(ii) = std(run.LeavingRR(run.PatchOrder == ii));
end

% how much reward did they get?
TotalReward = sum(run.Reward);

end


