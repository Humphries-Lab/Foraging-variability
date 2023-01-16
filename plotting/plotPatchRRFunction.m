%% plot reward rate function for patches

RR_Leave = [21.8678 18.5632]; % this is the average background RR at which subjects should leave, for Rich and Poor env... -
... respectively,  derived from running Nils' script with 32.5 45 57.5 as A, 0.075 as a, and a 6 sec delay between patches
    ... this value also corresponds to the instantaneous RR on the gain function at which subjects should leave.
    clear optLT
A=[32.5 45 57.5];a=0.075;
for e=1:2 %each env
    for p=1:3 % each patch type
        optLT(e,p)=(log(RR_Leave(e)/A(p)))/-a;
    end
end

environ.time_step = 1; % adjust time step if required

for j=1:length(A) % for each patch type
    fun = @(T) A(j)*exp(-T*0.075); %Exponential function currently
    rr(:,j) = fun(0:50);
    k = 0:0.05:50;
    for i = 1:length(k)
        r(i,j) = integral(fun,k(i),k(i)+environ.time_step); % This creates matrix 'q' with each cell reward earned in time_step epoch at that reward rate
    end
end

c = [0 0 0;0.4  0.4  0.4; 0.7    0.7    0.7];
figure
plot(0:50, rr(:,1), 'LineWidth', 2, 'Color', c(1,:))
hold on
plot(0:50, rr(:,2), 'LineWidth', 2, 'Color', c(2,:))
plot(0:50, rr(:,3), 'LineWidth', 2, 'Color', c(3,:))
scatter(optLT(2,:), RR_Leave(2), [], [0    0.4471    0.7412], 'filled')
scatter(optLT(1,:), RR_Leave(1), [], [0.6353    0.0784    0.1843], 'filled')
hold off

xlabel('Time in patch (s)', 'FontSize', 14)
ylabel('Foreground reward rate', 'FontSize', 14)
poor_ref = refline(0, environ.averageRR.poor);
rich_ref = refline(0, environ.averageRR.rich);

poor_ref.Color = [0    0.4471    0.7412];
rich_ref.Color = [0.6353    0.0784    0.1843];

poor_ref.LineStyle = '--';
rich_ref.LineStyle = '--';
title('Reward rate of different patch types', 'FontSize', 16)
legend('Low yield', 'Medium yield', 'High yield')

% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])
% print('/Users/exs165/Dropbox/foraging-project/results/rr_function', '-dsvg')