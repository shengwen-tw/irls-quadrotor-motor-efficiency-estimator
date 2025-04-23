% Noise levels
noise_levels = [0, 5, 10, 15, 20];

voltage_rmse = [
    0.000697, 0.000434, 0.000663, 0.000425; % 0% noise
    0.012480, 0.012171, 0.012358, 0.012283; % 5% noise
    0.026143, 0.025905, 0.025976, 0.025873; % 10% noise
    0.041932, 0.042572, 0.042239, 0.041778; % 15% noise
    0.058325, 0.055829, 0.056704, 0.057013  % 20% noise
    ];

fault_rmse = [
    0.006577, 0.006605, 0.006596, 0.006527; % 0% noise
    0.016067, 0.015575, 0.015908, 0.015716; % 5% noise
    0.028781, 0.029926, 0.031633, 0.028724; % 10% noise
    0.042765, 0.042744, 0.043298, 0.041626; % 15% noise
    0.058700, 0.056542, 0.054950, 0.055535  % 20% noise
    ];

combined_rmse = [
    0.016517, 0.015921, 0.016138, 0.016278; % 0% noise
    0.020052, 0.019914, 0.020053, 0.019977; % 5% noise
    0.029403, 0.027962, 0.028835, 0.029225; % 10% noise
    0.046326, 0.044790, 0.045545, 0.044808; % 15% noise
    0.058368, 0.057877, 0.056776, 0.057518  % 20% noise
    ];

% Define motor names
motor_names = arrayfun(@(i) sprintf('Motor %d', i), 1:4, 'UniformOutput', false);

% Create figure
fig = figure('Name', 'RMSE');
titles = {'Voltage degradation', 'Fault injection', 'Voltage degradation + Fault injection'};
datasets = {voltage_rmse, fault_rmse, combined_rmse};

for subplot_idx = 1:3
    subplot(3,1,subplot_idx);
    hold on;
    
    % Get current dataset
    data = datasets{subplot_idx};
    
    % Plot RMSE for each motor
    for motor_idx = 1:4
        plot(noise_levels, data(:, motor_idx), ...
            'LineWidth', 1.7, ...
            'Marker', 'd', ... % Diamond marker
            'DisplayName', motor_names{motor_idx});
    end
    
    % Compute dynamic y-ticks
    max_val = max(data(:));
    y_max = ceil(max_val * 100) / 100;        % Round up to nearest 0.01
    yticks(0:0.01:y_max);                     % Set fine ticks
    ylim([0 y_max]);                          % Match the y-axis limits
    
    title(titles{subplot_idx}, 'FontSize', 11);
    ylabel('RMSE', 'FontSize', 11);
    grid on;
    box on;
    
    lgd = legend('Location', 'northwest');
    lgd.NumColumns = 2;
    
    if subplot_idx == 3
        xlabel('Noise level [%]', 'FontSize', 11);
    end
end

sgtitle('Motor efficiency RMSE under different conditions');
exportgraphics(fig, 'RMSE.png');

disp("Press any key to leave");
pause;
close all;
