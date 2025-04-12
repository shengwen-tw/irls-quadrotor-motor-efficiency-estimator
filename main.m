math = se3_math;
estimator = motor_estimator;
rng(4419520);

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

x_avg = [1; 1; 1; 1];
alpha = 0.05;

time_arr = zeros(1, ITERATION_TIMES);
x_arr = zeros(4, ITERATION_TIMES);
x_error_arr = zeros(4, ITERATION_TIMES);

for i = 1:ITERATION_TIMES
    x = [1; 1; 1; 1];
    fprintf("iteration: %d\n", i)
    
    batch =  get_new_batch(data, i, estimator.n, 0);
    
    % Run motor efficiency estimator
    [x, skip] = estimator.run(i, batch, x);
    %if skip == 1
    %    continue;
    %end
    
    % Low-pass filter for x
    x_avg = alpha*x + (1-alpha)*x_avg;
    
    %disp(x_avg)
    disp(x)
    
    % Record
    time_arr(i) = i * dt;
    x_arr(:, i) = x;
    x_error_arr(:, i) = true_efficiency(:, i) - x;
end

% RMSE
rmse = sqrt(mean(x_error_arr.^2, 2));
fprintf("RMSE = (%f, %f, %f, %f)\n", rmse(1), rmse(2), rmse(3), rmse(4));

% Estimation result
figure('Name', 'Motor efficiency');
subplot (4, 1, 1);
plot(time_arr, x_arr(1, :));
hold on;
h2 = plot(time_arr, true_efficiency(1, 1:ITERATION_TIMES), 'r', 'LineWidth', 1.7);
legend(h2, 'true efficiency', 'Location', 'Best');
title('Motor efficiency');
xlabel('time [s]');
ylabel('\eta_1');
ylim([0.2 1.2]);
subplot (4, 1, 2);
plot(time_arr, x_arr(2, :))
hold on;
h2 = plot(time_arr, true_efficiency(2, 1:ITERATION_TIMES), 'r', 'LineWidth', 1.7);
legend(h2, 'true efficiency', 'Location', 'Best');
xlabel('time [s]');
ylabel('\eta_2');
ylim([0.2 1.2]);
subplot (4, 1, 3);
plot(time_arr, x_arr(3, :));
hold on;
h2 = plot(time_arr, true_efficiency(3, 1:ITERATION_TIMES), 'r', 'LineWidth', 1.7);
legend(h2, 'true efficiency', 'Location', 'Best');
xlabel('time [s]');
ylabel('\eta_3');
ylim([0.2 1.2]);
subplot (4, 1, 4);
plot(time_arr, x_arr(4, :));
hold on;
h2 = plot(time_arr, true_efficiency(4, 1:ITERATION_TIMES), 'r', 'LineWidth', 1.7);
legend(h2, 'true efficiency', 'Location', 'Best');
xlabel('time [s]');
ylabel('\eta_4');
ylim([0.2 1.2]);

% Error
figure('Name', 'Error of motor efficiency');
subplot (4, 1, 1);
plot(time_arr, x_error_arr(1, :));
title('Motor efficiency error');
xlabel('time [s]');
ylabel('\eta_1 error');
subplot (4, 1, 2);
plot(time_arr, x_error_arr(2, :));
xlabel('time [s]');
ylabel('\eta_2 error');
subplot (4, 1, 3);
plot(time_arr, x_error_arr(3, :));
xlabel('time [s]');
ylabel('\eta_3 error');
subplot (4, 1, 4);
plot(time_arr, x_error_arr(4, :));
xlabel('time [s]');
ylabel('\eta_4 error');

disp("Press any key to leave");
pause;
close all;

function batch = get_new_batch(data, iteration, len, random)
if random == 1
    idx = randi([1, ITERATION_TIMES]);
else
    idx = iteration;
end

idx_end = idx+len-1;
batch.p = data.pos_arr(:, idx:idx_end);
batch.v = data.vel_arr(:, idx:idx_end);
batch.W = data.W_arr(:, idx:idx_end);
batch.R = data.R_arr(:, :, idx:idx_end);
batch.f = data.f_arr(:, idx:idx_end);
batch.f_motors = data.f_motors_arr(:, idx_end);
batch.M = data.M_arr(:, idx:idx_end);
batch.m = data.m;
batch.J = data.J;
batch.c = data.c;
batch.d = data.d;
end