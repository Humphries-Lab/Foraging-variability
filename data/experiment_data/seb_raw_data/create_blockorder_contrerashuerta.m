clear
lt = readtable('../contrerashuerta_trialbytrial.csv','ReadVariableNames',true);

data = table;
for iS = [450:452,454:450+lt.sub(end)]
    tmp = load(sprintf('%d.mat',iS));
    subj_data=struct2table(tmp.result.data);
    subj_data = subj_data(subj_data.tEnd>0.5,:);

    subj_data = subj_data(~isoutlier(subj_data.tEnd,'mean'),:); %% don't use interquartlie range ('quartiles') as it's too strict. 
    subj_data = subj_data(subj_data.condN == 1,:); % restrict to self beneficiary only
    data = table(subj_data.block, subj_data.blVal);

    subjOrder = unique(data, 'rows','stable');
    blockOrder(iS-449,:) = subjOrder.Var2;
    switchIndex = diff(data.Var1) ~= 0;  
    blockSwitchIndex{iS-449} = [1; switchIndex]; % append with 1 (first patch in first block); 

end

blockOrder = blockOrder([1:3,5:end],:);
blockSwitchIndex = blockSwitchIndex([1:3,5:end]);

writematrix(blockOrder, 'contrerashuerta_blockOrder.csv')
save('contrerashuerta_blockSwitchIndex.mat','blockSwitchIndex')