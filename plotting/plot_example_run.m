function fig = plot_example_run(data, model)

if ismember(model, [1, 2])% if using MVT model

    %% calculate agent's averageRR based on their actions

    % if contains(block, 'rich')
    %     averageRR = 21.8678;
    % elseif contains(block, 'poor')
    %     averageRR = 18.5632;
    % end

    %% plot
    clf;
    fig = figure(1);
    fig.Position = [100 100 900 400];
    subplot(2,1,1)
    hold on
    plot(data.time_in_block, data.patchRR, 'LineWidth', 1) % plot actual patch rr over block time
    hline = refline(0,mean(data.patch_reward)); % plot average reward they obtained
    hline.Color = [0.8500 0.3250 0.0980];
    hline.LineWidth = 1;
    legend('Patch rr', 'Average rr', 'FontSize', 14, 'FontName', 'Helvetica', 'Location', 'southoutside')
    title(sprintf('M%d - Example simulation run', model),'FontSize', 18, 'FontName', 'Helvetica');
    xlabel('Time (s)','FontSize', 14, 'FontName', 'Helvetica');
    ylabel('Reward rate','FontSize', 14, 'FontName', 'Helvetica');

    subplot(2,1,2);
    hold on
    plot(data.time_in_block, data.patchRR, 'LineWidth', 1) % plot estimated patch rr over block time
    plot(data.time_in_block, data.rho, 'LineWidth', 1) % plot estimated average rr over block time
    legend('Estimated patch rr', '\rho: Estimated average rr', 'FontSize', 14, 'FontName', 'Helvetica', 'Location', 'southoutside')
    xlabel('Time (s)','FontSize', 14, 'FontName', 'Helvetica');
    ylabel('Reward rate','FontSize', 14, 'FontName', 'Helvetica');


else % if using any of the R-learning models

    %% plot
    clf;
    fig = figure(1);
    fig.Position = [100 100 1000 500];
    subplot(2,1,1);
    hold on
    plot(data.time_in_block, data.Q_leave, 'LineWidth', 1) % plot Q(leave) over block time
    plot(data.time_in_block, data.Q_stay, 'LineWidth', 1) % plot Q(stay) over block time
    legend('Q-leave', 'Q-stay', 'FontSize', 16, 'FontName', 'Helvetica')
    title(sprintf('M%d - Example simulation run', model),'FontSize', 18, 'FontName', 'Helvetica');
    xlabel('Time (s)','FontSize', 16, 'FontName', 'Helvetica');
    ylabel('Q','FontSize', 16, 'FontName', 'Helvetica');

    subplot(2,1,2)
    hold on
    plot(data.time_in_block, data.patchRR, 'LineWidth', 1) % plot patch reward rate
    plot(data.time_in_block, data.rho*57.5, 'LineWidth', 1) % plot estimated average RR
    legend('Patch RR', '\rho: Estimated average RR', 'FontSize', 16, 'FontName', 'Helvetica')
    xlabel('Time (s)','FontSize', 16, 'FontName', 'Helvetica');
    ylabel('Reward rates','FontSize', 16, 'FontName', 'Helvetica');

end

end



