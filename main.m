math = se3_math;
estimator = motor_estimator;
rng(4419520);

data = load('sim_log.mat');
data_size = size(data.time_arr, 2);
fprintf('Data size = %d\n', data_size);

dt = data.dt;
g = data.g;
m = data.m;
J = data.J;
e3 = [0; 0; 1];
I_3x3 = eye(3, 3);

x_pred = zeros(3, data_size);
v_pred = zeros(3, data_size);
W_pred = zeros(3, data_size);
R_pred = zeros(3, 3, data_size);
euler_pred = zeros(3, data_size);

% Initialize the first step
x_pred(:, 1) = data.pos_arr(:, 1);
v_pred(:, 1) = data.vel_arr(:, 1);
W_pred(:, 1) = data.W_arr(:, 1);
R_pred(:, :, 1) = data.R_arr(:, :, 1);
euler_pred(:, 1) = data.euler_arr(:, 1);

%%%%%%%%%%%%%%%%%%%%%%
% Run ADMM estimator %
%%%%%%%%%%%%%%%%%%%%%%
fprintf("Before estimation\n");
fprintf("m = %f\n", data.m);
fprintf("Jx = %f\n", data.J(1, 1));
fprintf("Jy = %f\n", data.J(2, 2));
fprintf("Jz = %f\n", data.J(3, 3));

estimator = estimator.init(math, dt);
x_est = estimator.run(data);

disp("Press any key to leave");
pause;
close all;