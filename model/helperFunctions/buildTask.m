function task = buildTask(study)

% Foraging variability project 
% Emma Scholey, 16 January 2025 

% BUILDTASK generates task variables for each dataset

% INPUTS:
%   study: either field-human, berry-human, or rat 

switch study

    case 'leheron'

        % task parameters
        task.travelTime = [6, 6]; % delay between patches. same time in both environments
        task.timeStep = 1;  % T, seconds (discretising continuous time)
        task.blockTime = 600; % seconds in block
        task.r0 = [32.5, 45, 57.5]; % starting reward rate
        task.patchNames = {'Low', 'Mid', 'High'};

        task.decayRate = 0.075; % decay rate
        task.rewardFunction = 'exponential';
        
        task.environNames = {'rich', 'poor'};
        task.nEnviron = 2;
        task.nBlocks = 2;
        task.nPatch = numel(task.r0);

        % specify which patch the agent enters. This trial order was set by
        % CLH for cabergoline data, where every 10 patches leads to correct
        % proportions of low, medium or high patches based on the environment type

        % patch order in rich environment
        task.patchOrder{1} = repmat([3,2,3,1,3,2,1,3,3,2,1,3,2,1,3,3,2,3,3,2,3,3,1,3,3,2,3,1,2,2,3,2,2,1,3,3,2,3,1,3,3,3,2,2,3,3,3,1,1,2,3,2,2,3,3,1,1,3,2,3,3,2,3,1,1,3,2,2,3,3,3,3,3,1,3,2,3,2,2,1,1,3,2,3,3,1,2,2,3,3,2,3,1,3,3,3,1,2,2,3], [1,5]);
        % patch order in poor environment
        task.patchOrder{2} = repmat([1,1,2,1,3,2,2,1,3,1,1,2,3,1,2,1,1,2,3,1,1,1,2,3,1,2,3,2,1,1,3,1,1,2,2,1,2,3,1,1,1,2,3,1,1,2,1,3,1,2,3,2,2,1,2,1,1,1,3,1,1,1,2,1,2,1,1,2,3,3,1,3,1,3,2,2,2,1,1,1,2,3,1,3,1,2,1,2,1,1,2,1,2,3,2,3,1,1,1,1], [1,5]);

        task.optAvgRR(1) = simulateOptimalRho(task.decayRate, task.rewardFunction, task.r0, task.patchOrder{1},task.travelTime(1), task.blockTime); % 21.25 
        task.optAvgRR(2) = simulateOptimalRho(task.decayRate, task.rewardFunction, task.r0, task.patchOrder{2},task.travelTime(2), task.blockTime); % 18.30

        % convert reward rates into optimal leaving times as defined by mvt
        for e=1:task.nEnviron %each env
            for p=1:task.nPatch % each patch type
                task.optLT(e,p)=(log(task.optAvgRR(e)/task.r0(p)))/-task.decayRate;
            end
        end

    case 'contrerashuerta'

        % task parameters
        task.travelTime = [3, 5]; % delay in s between patches [rich env, poor env] 
        task.timeStep = 1;  % T, seconds (discretising continuous time)
        task.blockTime = 300; % seconds in block
        task.r0 = [34.5, 57.5]; % starting reward rate
        task.patchNames = {'Low', 'Mid'};

        task.decayRate = 0.11; % decay rate
        task.rewardFunction = 'exponential';

        task.environNames = {'rich', 'poor'};
        task.nEnviron = 2;
        task.nBlocks = 6;
        task.nPatch = numel(task.r0);

        % patch order in rich environment - note, provide more patches than
        % we actually require. since don't know how many patches the agent
        % will visit in their alloted time! 
        task.patchOrder{1} = repmat([2,2,1,1,2,1,2,1,2,1], [1,40]); % assume same random patch order in both environments (since environment determined by travel time, not patch likelihood)
        % patch order in poor environment
        task.patchOrder{2} = repmat([2,2,1,1,2,1,2,1,2,1], [1,40]);

        task.optAvgRR(1) = simulateOptimalRho(task.decayRate, task.rewardFunction, task.r0, task.patchOrder{1},task.travelTime(1), task.blockTime); % 23.09
        task.optAvgRR(2) = simulateOptimalRho(task.decayRate, task.rewardFunction, task.r0, task.patchOrder{2},task.travelTime(2), task.blockTime); % 19.10

        % convert reward rates into optimal leaving times as defined by mvt
        for e=1:task.nEnviron %each env
            for p=1:task.nPatch % each patch type
                task.optLT(e,p)=(log(task.optAvgRR(e)/task.r0(p)))/-task.decayRate;
            end
        end

    case 'kane'
         % task parameters. Note that 1 state = 10s (trial by trial rather than
         % continuous time) 
        task.travelTime = [1, 3]; % delay in STATES between patches [rich env, poor env] 
        task.timeStep = 1;  % T, seconds (discretising continuous time). Here this is equivalent to one state
        task.blockTime = 300; % how many trials in each environment (for simulation only) 
        task.r0 = [60, 90, 120]; % starting reward rate
        task.patchNames = {'Low', 'Mid', 'High'};

        task.decayRate = 8; % decay rate
        task.rewardFunction = 'linear'; % linear or exponential

        task.environNames = {'rich', 'poor'};
        task.nEnviron = 2;
        task.nBlocks = 10;
        task.nPatch = numel(task.r0);

        % patch order in rich environment
        task.patchOrder{1} = repmat([2 3 1 2 1 3 1 2 3 1 3 2 1 2 3 3 1 2 1 3 2], [1,40]); % assume same random patch order in both environments
        % patch order in poor environment
        task.patchOrder{2} = repmat([2 3 1 2 1 3 1 2 3 1 3 2 1 2 3 3 1 2 1 3 2], [1,40]);

        task.optAvgRR(1) = simulateOptimalRho(task.decayRate, task.rewardFunction, task.r0, task.patchOrder{1},task.travelTime(1), task.blockTime); % 66.5
        task.optAvgRR(2) = simulateOptimalRho(task.decayRate, task.rewardFunction, task.r0, task.patchOrder{2},task.travelTime(2), task.blockTime); % 46.49

        % convert reward rates into optimal leaving times as defined by mvt
        for e=1:task.nEnviron %each env
            for p=1:task.nPatch % each patch type
                task.optLT(e,p)=(task.optAvgRR(e)-task.r0(p))/-task.decayRate;
            end
        end

        task.optLT(1,1) = 1; % replace optimal time for harvests for richest patch in richest environment to 1s - mouse needs to make one leverpress first 
       
        % since showing harvests per timestep, round optimal leave decision
        % up to nearest time step 
        task.optLT = round(task.optLT);

end
