
clear
close all

%% Read model and simulation tables

File = '/Users/exs165/Dropbox/foraging/code/Francisca/230720_FranciscaDevelopment/Model_options.xlsx';

Model_table = readtable(File,'Sheet','Model');
SimTable= readtable(File,'Sheet','Simulations');

% Load Simulations to run
RunSim = Load_SimulationParameters(SimTable, Model_table);

% Clear the tables not needed anymore
clear SimTable Model_table

%% Run Simulations

Results= cell([height(RunSim) 1]);
LT=cell([height(RunSim) 1]);

for i = 1:height(RunSim)
    fprintf('Simulation progress: %d/%d\n',i,height(RunSim));
    input = RunSim(i,:);    
    out = Simulate_RLOn_Model(input.RLparams,input.BlockOrder,input.Task,input.Model);
           
    % Extract leaving times (LT) for each patch and block type 
    for patch = 1: RunSim.Task(i).PatchNum
        for block = 1: RunSim.Task(i).BlockNum
            LT{i}(patch,block) = mean(out.LeavingTime{block}(out.PatchOrder{block}(1:end-1) == patch), 'omitnan');
        end
    end
    % Store results
    Results{i} = out; 
    
end

%% Store results in last columns of Simulation table   
RunSim = addvars(RunSim,Results,LT,'after','Model');

%Clear vbs
clear Results LT

%% compute means and standard errors across simulations IDs

%Group simulations by SimID
ID = unique(RunSim.SimID);
%unique([RunSim.Model(1:4).Num]);

% Initialise vbs to store results
LTMean= cell([length(ID) 1]);
LTSEM= cell([length(ID) 1]);

for n=1:length(ID)
    %find number and position of simulations with same ID
    pos= find(RunSim.SimID == ID(n));   
    NSim =sum(RunSim.SimID == ID(n));
    % compute means and standard errors across simulations
    tempM= reshape(horzcat(RunSim.LT{pos(1):pos(end)}),3,2,[]);
    LTMean{n} = squeeze(mean(tempM,3));
    LTSEM{n} = squeeze(std(tempM,[],3))./sqrt(NSim);     
end

clear tempM NSim pos n

%% Plot mean leaving times

close all

% set up colours
c = [0.6353    0.0784    0.1843; 0    0.4471    0.7412]; % [rich, poor]
% %labels
% xlabel('Patch quality')
% ylabel('Patch leaving times (s)')

for n=1:length(ID)
    
    %Create new figure
    figure(n)
    hold on
    %Student's t inverse cumulative distribution function
    NSim= sum(RunSim.SimID == ID(n));
    ts = tinv([0.025  0.975],NSim-1);      % T-Score
    LTCI = LTSEM{n}.* ts(2);    
    %line plot with error bars
    errorbar([1 2 3], LTMean{n}(:,1), LTCI(:,1), 'LineWidth', 2, 'Color', c(1,:)); % rich environment
    errorbar([1 2 3], LTMean{n}(:,2), LTCI(:,2), 'LineWidth', 2, 'Color', c(2,:)); % poor environment
    
    % tidy up axes
    pos= find(RunSim.SimID == ID(n));
    title(sprintf('Model n.%d', RunSim.Model(pos(1)).Num))
    ax = gca;
    xlim(ax, xlim(ax) + [-1,1]*range(xlim(ax)).* 0.1)
    xticks(ax, [1 2 3]); xticklabels(ax, {'Low'; 'Medium'; 'High'}); xlabel('Patch quality');
    ylim([0, 30]); ylabel('Patch leaving times (s)');
    
    %legend
    %legend({'Subject - rich', 'Subject - poor', 'Agent - rich', 'Agent - poor'}, 'Location', 'north');
    legend({'Agent - rich', 'Agent - poor'}, 'Location', 'north');
    %subtitle(sprintf('Error bars show 95%% CI. NSim = %d', NSim))

end













