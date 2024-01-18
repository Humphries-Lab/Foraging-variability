clear
lt = readtable('leheron_trialbytrial.csv','ReadVariableNames',true);

for iS = 1:39

    their_data = lt(lt.sub==iS,:);

    switchIndex = diff(their_data.env) ~= 0;
    blockSwitchIndex{iS} = [1; switchIndex]; % append with 1 (first patch in first block);
end

save('leheron_blockSwitchIndex.mat','blockSwitchIndex')