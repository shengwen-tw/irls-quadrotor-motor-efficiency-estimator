math = se3_math;
estimator = motor_estimator;

data = load('sim_log.mat');
data_size = size(data.time_arr, 2);
fprintf('Data size = %d\n', data_size);

true_efficiency = data.motor_efficiency_arr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run motor efficiency estimation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = data.dt;
ITERATION_TIMES = data_size - estimator.n + 1;
estimator = estimator.init(math, dt);

time_arr = zeros(1, ITERATION_TIMES);
x_arr = zeros(4, ITERATION_TIMES);
x_error_arr = zeros(4, ITERATION_TIMES);
cond_arr = zeros(1, ITERATION_TIMES);

x_last = [1; 1; 1; 1];
for i = 1:ITERATION_TIMES
    fprintf("iteration: %d\n", i)
    
    batch =  get_new_batch(data, i, estimator.n);
    
    % Run motor efficiency estimator
    [x, cond] = estimator.run(i, batch, [0.5; 0.5; 0.5; 0.5], x_last);
    x_last = x;
    
    disp(x)
    
    % Record
    time_arr(i) = i * dt;
    x_arr(:, i) = x;
    x_error_arr(:, i) = true_efficiency(:, i) - x;
    cond_arr(i) = cond;
end

% RMSE / Varaince / Standard deviation
rmse = sqrt(mean(x_error_arr.^2, 2));
variance = var(x_error_arr, 0, 2);
std_dev = std(x_error_arr, 0, 2);
fprintf("RMSE = (%f, %f, %f, %f)\n", rmse(1), rmse(2), rmse(3), rmse(4));
fprintf("Variance = (%f, %f, %f, %f)\n", variance(1), variance(2), variance(3), variance(4));
fprintf("Standard Deviation = (%f, %f, %f, %f)\n", std_dev(1), std_dev(2), std_dev(3), std_dev(4));

% Estimation result
figure('Name', 'Motor efficiency: Estimated vs True');
subplot (4, 1, 1);
h1 = plot(time_arr, x_arr(1, :), 'LineWidth', 1.7);
hold on;
h2 = plot(time_arr, true_efficiency(1, 1:ITERATION_TIMES), 'r', 'LineWidth', 1.7);
legend('Estimated', 'True', 'Location', 'Best', 'Orientation', 'horizontal');
title('Motor efficiency: Estimated vs True');
xlabel('time [s]');
ylabel('\eta_1');
ylim([0.2 1.2]);
grid on;
subplot (4, 1, 2);
h1 = plot(time_arr, x_arr(2, :), 'LineWidth', 1.7);
hold on;
h2 = plot(time_arr, true_efficiency(2, 1:ITERATION_TIMES), 'r', 'LineWidth', 1.7);
legend('Estimated', 'True', 'Location', 'Best', 'Orientation', 'horizontal');
xlabel('time [s]');
ylabel('\eta_2');
ylim([0.2 1.2]);
grid on;
subplot (4, 1, 3);
h1 = plot(time_arr, x_arr(3, :), 'LineWidth', 1.7);
hold on;
h2 = plot(time_arr, true_efficiency(3, 1:ITERATION_TIMES), 'r', 'LineWidth', 1.7);
legend('Estimated', 'True', 'Location', 'Best', 'Orientation', 'horizontal');
xlabel('time [s]');
ylabel('\eta_3');
ylim([0.2 1.2]);
grid on;
subplot (4, 1, 4);
h1 = plot(time_arr, x_arr(4, :), 'LineWidth', 1.7);
hold on;
h2 = plot(time_arr, true_efficiency(4, 1:ITERATION_TIMES), 'r', 'LineWidth', 1.7);
legend('Estimated', 'True', 'Location', 'Best', 'Orientation', 'horizontal');
xlabel('time [s]');
ylabel('\eta_4');
ylim([0.2 1.2]);
grid on;

% Error
figure('Name', 'Error of motor efficiency');
subplot (4, 1, 1);
plot(time_arr, x_error_arr(1, :), 'LineWidth', 1.7);
title('Motor efficiency error');
xlabel('time [s]');
ylabel('\eta_1 error');
grid on;
subplot (4, 1, 2);
plot(time_arr, x_error_arr(2, :), 'LineWidth', 1.7);
xlabel('time [s]');
ylabel('\eta_2 error');
grid on;
subplot (4, 1, 3);
plot(time_arr, x_error_arr(3, :), 'LineWidth', 1.7);
xlabel('time [s]');
ylabel('\eta_3 error');
grid on;
subplot (4, 1, 4);
plot(time_arr, x_error_arr(4, :), 'LineWidth', 1.7);
xlabel('time [s]');
ylabel('\eta_4 error');
grid on;

figure('Name', 'Condition number over time');
%
y_thresh = estimator.cond_threshold;
x_fill = [time_arr(1), time_arr(end), time_arr(end), time_arr(1)];
y_fill = [y_thresh, y_thresh, max(cond_arr)*1.1, max(cond_arr)*1.1];
fill(x_fill, y_fill, [1 0.8 0.8], 'EdgeColor', 'none');  % Light red fill
hold on;
%
%plot(time_arr, cond_arr, 'LineWidth', 1.7, 'Color', '#1f77b4');
bar(time_arr, cond_arr, 'FaceColor', '#1f77b4');
h2 = yline(estimator.cond_threshold, '--k', 'Color', 'r', 'LineWidth', 1.7);
legend(h2, 'Rejecting threshold');
xlabel('time [s]');
ylabel('Condition number of J_f^TGJ_f');
xlim([time_arr(1), time_arr(end)]);
%ylim([0, 1/estimator.cond_threshold * 1.01])
set(gca, 'YScale', 'log')
grid on;

disp("Press any key to leave");
pause;
close all;

function batch = get_new_batch(data, iteration, len)
idx = iteration;

idx_end = idx+len-1;
batch.p = data.pos_arr(:, idx:idx_end);
batch.v = data.vel_arr(:, idx:idx_end);
batch.W = data.W_arr(:, idx:idx_end);
batch.R = data.R_arr(:, :, idx:idx_end);
batch.f = data.f_arr(:, idx:idx_end);
batch.f_motors = data.f_motors_arr(:, idx:idx_end);
batch.motor_efficiency = data.motor_efficiency_arr(:, idx_end); % Ground truth
batch.M = data.M_arr(:, idx:idx_end);
batch.m = data.m;
batch.J = data.J;
batch.c = data.c;
batch.d = data.d;
end