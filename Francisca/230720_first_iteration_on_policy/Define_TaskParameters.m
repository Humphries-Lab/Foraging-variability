
function [Task] = Define_TaskParameters(Study, BlockPresentation)

switch Study
    case 'LeHeron'
       
        % TASK STRUCTURE
        Task.TimeStep=1;     % time in seconds to update time_in_block and time_in_patch
        Task.TravelTime = 6; % delay between patches in seconds
        Task.BlockNum=2;     % Number of blocks in the task (Block=1 Rich Env.; Block=2 Poor Env)
        Task.BlockTime=600;  % Duration of each block in seconds

        switch BlockPresentation    % Total duration of the task in seconds depending on block presentation
           case 'Combined' 
           Task.Duration= Task.BlockTime* Task.BlockNum;    
            otherwise
             Task.Duration= Task.BlockTime;      
        end
         
        % PATCH TYPES
        Task.PatchNum =3;                           % Number of patch types (integer [1,2,..,N])
        InitialYield = [32.5, 45, 57.5];            % Initial Patch Yield (S) in LeHeron task
        t= 0:Task.TimeStep:Task.Duration;           % time vector until end of the task
        for p=1:Task.PatchNum                       % for each yield type
            RRfun = @(T) (InitialYield(p)*exp(-T*0.075));                        % Patch Foreground reward rate function
            %Patch.RR(:,p) = RRfun(0:TaskTime);
            for epoch=1:length(t)
                Task.PatchReward(epoch,p)= integral(RRfun,t(epoch),t(epoch)+Task.TimeStep); % Observed reward per time_step epoch at patch RR            
            end
        end
        
        % BACKGROUND REWARD RATE PER ENVIRONMENT (PATCH TYPE PROPORTION)
        % Patch Order in Rich Environment
        Task.Block_POrder(:,1) = repmat([3,2,3,1,3,2,1,3,3,2,1,3,2,1,3,3,2,3,3,2,3,3,1,3,3,2,3,1,2,2,3,2,2,1,3,3,2,3,1,3,3,3,2,2,3,3,3,1,1,2,3,2,2,3,3,1,1,3,2,3,3,2,3,1,1,3,2,2,3,3,3,3,3,1,3,2,3,2,2,1,1,3,2,3,3,1,2,2,3,3,2,3,1,3,3,3,1,2,2,3], [1,5]);
        % Patch Order in Poor Environment
        Task.Block_POrder(:,2) = repmat([1,1,2,1,3,2,2,1,3,1,1,2,3,1,2,1,1,2,3,1,1,1,2,3,1,2,3,2,1,1,3,1,1,2,2,1,2,3,1,1,1,2,3,1,1,2,1,3,1,2,3,2,2,1,2,1,1,1,3,1,1,1,2,1,2,1,1,2,3,3,1,3,1,3,2,2,2,1,1,1,2,3,1,3,1,2,1,2,1,1,2,1,2,3,2,3,1,1,1,1], [1,5]);

        % MVT PREDICTION
        Task.MTVOptimal(1) = 21.8678; % optimal average RR in rich environment
        Task.MTVOptimal(2) = 18.5632; % optimal average RR in poor environment        
        
        
    otherwise        
        % Specify here other studies' data or custom parameters for the
        % task:
        % TimeStep, BlockNum, BlockTime, Duration, PatchNum,PatchReward,Block_POrder 
           
        
end
end