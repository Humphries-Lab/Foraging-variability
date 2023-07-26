
function [Simulations] = Load_SimulationParameters(SimTable, ModelTable)

% This function creates a Simulation table that combines the infomation from the excel file (it indicates the model type, the
% study , and the number of simulations to run).
% 
%   INPUTS:
%   SimTable: A table that defines the simulations to run and contains the following information: 
%           -Study to simulate which will define the task parameters
%           -Number of Model to run (according to ModelTable)
%           -RL parameters: AlphaQ, AlphaRho, Beta, and Epsilon (in case of e-greedy model)
%           -Number of Simulations to run with the specified model
%           -BlockOrder: the order in which the different types of
%           environments are presented (options are: alternate each order for 50% of simulations, Rich-poor, Poor-rich, random order)
%           -Block Presentation: whether blocks/environments are presented
%           separately or in one combined session
%
%   ModelTable: A table that defines the model to run and has these parameters:
%           -Model number: It identifies the model
%           -Patch_sequence: Indicates whether Qstay sequence is multiple(per patch type) or single (regardless of patch type)
%           -Leave_sequence: Indicates whether Qleave has states (one per travel time) or is stateless (single value updated in each travel time)
%           -Next_patch: Prediccion of next patch (only if multiple patch
%           sequence). Types of prediction considered are 'previous'(same as patch just left), 
%           'Godlike'(next patch is known),'Stochastic'(guess randomly or based on the most frequent patch encountered in the block)
%           -Action_selection: Indicates whether next action is selected based on e-greedy or softmax (stochastic)function.
%
%   OUTPUT:
%       Simulations:   A table with as many rows as simulations to run that contains the following information:
%           -Study that indicates the task to simulate.
%           -RL parameters [Alpha Q, AlphaRho, Beta or Epsilon]
%           -BlockOrder: the order in which the environments will be presented
%           -Task parameters
%           -Model options
%
%   Version: 29/06/2023
%   Francisca Perea Perez


% READ TABLES FROM FILE IF THEY ARE NOT PROVIDED AS INPUT
if nargin < 2
    Tablefile = 'C:\Users\lpzfp\OneDrive - The University of Nottingham\FPP Documents\Matlab code\My_code\Model_options.xlsx';
    ModelTable = readtable(Tablefile,'Sheet','Model');
elseif nargin < 1
    Tablefile = 'C:\Users\lpzfp\OneDrive - The University of Nottingham\FPP Documents\Matlab code\My_code\Model_options.xlsx';
    SimTable= readtable(Tablefile,'Sheet','Simulations');
    ModelTable = readtable(Tablefile,'Sheet','Model');
end

% SET UP STORAGE FOR OUTPUT SIMULATION TABLE
Simulations= table(); 
Study =string.empty;
RLparam = double.empty;
BlockOrder=[];


% CREATE SIMULATION TABLE
for i=1: height(SimTable)       % loop across the simulations specified (row) in SimTable

   %First column:Simulation ID
   SimID= SimTable.SimID(i);
   
   %Second colum: Study
   Study = SimTable.Study{i};

   %Third column: RL parameters
   % Find the position in ModelTable of the model that we want to use in the
   % simulation
   pos= find(ModelTable.Model== SimTable.Model(i));  %find position in model table

   switch ModelTable.Action_selection{pos}          %if Stochastic selection I need Beta, if e-greedy I need Epsilon.
    case 'Stochastic'
        RLparam= [SimTable.AlphaQ(i) SimTable.AlphaRho(i) SimTable.Beta(i)];  % [AlphaQ AlphaRho Beta]
    otherwise
        RLparam= [SimTable.AlphaQ(i) SimTable.AlphaRho(i) (SimTable.Epsilon(i))/100]; % [AlphaQ AlphaRho Epsilon] Epsilon is introduced as percent probability
    end
    
   % Fourth column: Block order (the order in which both environments are presented)
    switch SimTable.BlockOrder{i}
        case 'Fifty_Fifty'          % 50% of the simulations per each order
            BlockOrder = [repmat([1 2], [SimTable.Nsimulations(i)/2, 1]); repmat([2 1], [SimTable.Nsimulations(i)/2, 1])];
        case 'Rich_Poor'            % Rich environment presented first, then poor 
            BlockOrder = repmat([1 2], [SimTable.Nsimulations(i), 1]);
        case 'Poor_Rich'            % Poor environment presented first, then rich
            BlockOrder = repmat([2 1], [SimTable.Nsimulations(i), 1]);        
        otherwise
            random = discreteinvrnd([0.5 0.5], SimTable.Nsimulations(i), 1);   % random order    
            BlockOrder(random==1,:)= repmat([1 2],length(random(random==1)),1);
            BlockOrder(random==2,:)= repmat([2 1],length(random(random==2)),1);
    end

    % Fifth column: Task parameters
    % Load parameters from the study/task we want to simulate
    Task = Define_TaskParameters(SimTable.Study{i}, SimTable.Block_presentation{i});

    % Sixth column: Model type
    % Store the characteristics of the model we want to simulate
    Model.Num= SimTable.Model(i);                               % Model Number 
    Model.Patch_sequence= ModelTable.Patch_sequence{pos};       % Patch sequence: either one per patch type or a single Qstay sequence for all.
    Model.Leave_sequence= ModelTable.Leave_sequence{pos};       % Leave sequence: either one state per travel time or a single stay for all.
    Model.Next_patch= ModelTable.Next_patch{pos};               % Prediction of next patch (only if multiple patch sequence).
    Model.Memory= ModelTable.Memory_NextPatch{pos};             % The memory of previous patches used to predict next patch (only if prediction is stochastic or egreedy)
    Model.Action_selection= ModelTable.Action_selection{pos};   % Action selection: either e-greedy or stochastic
    Model.Block_presentation= SimTable.Block_presentation{i};   % Simulate blocks/environments separately or in one combined session

    % Create temporal table with all parameters of this simulation
    tempTable= table(repmat(SimID,[SimTable.Nsimulations(i),1]),repmat(Study,[SimTable.Nsimulations(i),1]), repmat(RLparam,[SimTable.Nsimulations(i),1]),BlockOrder,repmat(Task,[SimTable.Nsimulations(i),1]), repmat(Model,[SimTable.Nsimulations(i),1]),...
        'VariableNames',{'SimID','Study','RLparams','BlockOrder','Task','Model'});

    % Add this simulation to the Simulations table
    Simulations= [Simulations; tempTable];

    % Clear variables
    tempTable=table();
    Study =string.empty;
    RLparam = double.empty;
    BlockOrder=[];    

end

end