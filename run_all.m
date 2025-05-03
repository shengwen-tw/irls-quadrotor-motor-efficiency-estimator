math = se3_math;

% Run `quadrotor_sim.m` to generate dataset first

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EKF: Baseline, has larger spikes                %
% Non-linear Least Square method: Has less spikes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = load('sim_log.mat');
data_size = size(data.time_arr, 2);
fprintf('Data size = %d\n', data_size);

true_efficiency = data.motor_efficiency_arr;  % 4 x T

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run motor efficiency estimation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = data.dt;

%======================
% Prepare two estimators
%======================
ekf = motor_estimator_ekf;
ekf = ekf.init(math, dt);
ITER_EKF = data_size - ekf.n + 1;

irls = motor_estimator_irls;
irls = irls.init(math, dt);
ITER_IRLS = data_size - irls.n + 1;

%======================
% Buffers: EKF
%======================
time_arr_ekf = (1:ITER_EKF) * dt;
x_arr_ekf = zeros(4, ITER_EKF);
x_err_ekf = zeros(4, ITER_EKF);
cond_ekf  = zeros(1, ITER_EKF);

%======================
% Buffers: IRLS
%======================
time_arr_irls = (1:ITER_IRLS) * dt;
x_arr_irls = zeros(4, ITER_IRLS);
x_err_irls = zeros(4, ITER_IRLS);
cond_irls  = zeros(1, ITER_IRLS);

%======================
% Run EKF loop
%======================
x_last = [1;1;1;1];
fprintf("\n=== Running EKF ===\n");
for i = 1:ITER_EKF
    fprintf("EKF iteration: %d/%d\n", i, ITER_EKF);
    batch = get_new_batch(data, i, ekf.n);
    [x, condv] = ekf.run(batch, x_last, x_last);
    x_last = x;
    
    x_arr_ekf(:, i) = x;
    x_err_ekf(:, i) = true_efficiency(:, i) - x;
    cond_ekf(i)     = condv;
end

%======================
% Run IRLS loop
%======================
x_last = [1;1;1;1];
fprintf("\n=== Running IRLS ===\n");
for i = 1:ITER_IRLS
    fprintf("IRLS iteration: %d/%d\n", i, ITER_IRLS);
    batch = get_new_batch(data, i, irls.n);
    [x, condv] = irls.run(batch, x_last, x_last);
    x_last = x;
    
    x_arr_irls(:, i) = x;
    x_err_irls(:, i) = true_efficiency(:, i) - x;
    cond_irls(i)     = condv;
end

%======================
% Metrics (EKF / IRLS)
%======================
rmse_ekf = sqrt(mean(x_err_ekf.^2, 2));
var_ekf  = var(x_err_ekf, 0, 2);
std_ekf  = std(x_err_ekf, 0, 2);

rmse_irls = sqrt(mean(x_err_irls.^2, 2));
var_irls  = var(x_err_irls, 0, 2);
std_irls  = std(x_err_irls, 0, 2);

%===========================================%
% Plots: two columns (left=EKF, right=IRLS) %
%===========================================%

% 1) Estimation result: 4 rows x 2 cols
fig_eff = figure('Name', 'Motor efficiency: Estimated vs True (EKF vs IRLS)');
tiledlayout(4,2, 'Padding','compact', 'TileSpacing','compact');

for m = 1:4
    % Left: EKF
    nexttile;
    plot(time_arr_ekf, x_arr_ekf(m, :), 'LineWidth', 1.7); hold on;
    plot(time_arr_ekf, true_efficiency(m, 1:ITER_EKF), 'r', 'LineWidth', 1.7);
    if m==1
        title('EKF: Estimated vs True', 'FontSize', 11);
        legend('Estimated','True', 'Location','Best', 'Orientation','horizontal');
    end
    ylabel(sprintf('\\eta_%d', m), 'FontSize', 11);
    if m==4, xlabel('time [s]', 'FontSize', 11); end
    ylim([0.2 1.2]); grid on;
    
    % Right: IRLS
    nexttile;
    plot(time_arr_irls, x_arr_irls(m, :), 'LineWidth', 1.7); hold on;
    plot(time_arr_irls, true_efficiency(m, 1:ITER_IRLS), 'r', 'LineWidth', 1.7);
    if m==1
        title('IRLS: Estimated vs True', 'FontSize', 11);
        legend('Estimated','True', 'Location','Best', 'Orientation','horizontal');
    end
    if m==4, xlabel('time [s]', 'FontSize', 11); end
    ylim([0.2 1.2]); grid on;
end
exportgraphics(fig_eff, 'estimation_result_ekf_irls.png');

% 2) Error: 4 rows x 2 cols
fig_err = figure('Name', 'Error of motor efficiency (EKF vs IRLS)');
tiledlayout(4,2, 'Padding','compact', 'TileSpacing','compact');

for m = 1:4
    % Left: EKF
    nexttile;
    plot(time_arr_ekf, x_err_ekf(m, :), 'LineWidth', 1.7);
    if m==1, title('EKF error'); end
    ylabel(sprintf('\\eta_%d error', m));
    if m==4, xlabel('time [s]'); end
    grid on;
    
    % Right: IRLS
    nexttile;
    plot(time_arr_irls, x_err_irls(m, :), 'LineWidth', 1.7);
    if m==1, title('IRLS error'); end
    if m==4, xlabel('time [s]'); end
    grid on;
end
exportgraphics(fig_err, 'error_ekf_irls.png');

% 3) Condition number: 1 row x 2 cols
fig_cond = figure('Name', 'Condition number over time (EKF vs IRLS)');
tiledlayout(1,2, 'Padding','compact', 'TileSpacing','compact');

% Left: EKF
nexttile;
y_thresh = ekf.cond_threshold;
x_fill = [time_arr_ekf(1), time_arr_ekf(end), time_arr_ekf(end), time_arr_ekf(1)];
y_fill = [y_thresh, y_thresh, max(cond_ekf)*1.1, max(cond_ekf)*1.1];
fill(x_fill, y_fill, [1 0.8 0.8], 'EdgeColor','none'); hold on;
bar(time_arr_ekf, cond_ekf, 'FaceColor', '#1f77b4');
h2 = yline(ekf.cond_threshold, '--k', 'Color', 'r', 'LineWidth', 1.7);
legend(h2, 'Rejecting threshold');
xlabel('time [s]', 'FontSize', 11);
ylabel('Condition number of J_f^TGJ_f', 'FontSize', 11);
xlim([time_arr_ekf(1), time_arr_ekf(end)]);
set(gca, 'YScale', 'log'); grid on;
title('EKF');

% Right: IRLS
nexttile;
y_thresh = irls.cond_threshold;
x_fill = [time_arr_irls(1), time_arr_irls(end), time_arr_irls(end), time_arr_irls(1)];
y_fill = [y_thresh, y_thresh, max(cond_irls)*1.1, max(cond_irls)*1.1];
fill(x_fill, y_fill, [1 0.8 0.8], 'EdgeColor','none'); hold on;
bar(time_arr_irls, cond_irls, 'FaceColor', '#1f77b4');
h2 = yline(irls.cond_threshold, '--k', 'Color', 'r', 'LineWidth', 1.7);
legend(h2, 'Rejecting threshold');
xlabel('time [s]', 'FontSize', 11);
ylabel('Condition number of J_f^TGJ_f', 'FontSize', 11);
xlim([time_arr_irls(1), time_arr_irls(end)]);
set(gca, 'YScale', 'log'); grid on;
title('IRLS');

exportgraphics(fig_cond, 'cond_ekf_irls.png');

%=========================================================================%
% Spikes in others-only fault windows (±3% band + time padding)           %
% Spikes = max(|x_err_*|) measured ONLY in others' windows (exclude self) %
%=========================================================================%

% Align to common length
T_common = min(size(x_err_ekf,2), size(x_err_irls,2));
T        = T_common;

% Inputs (computed earlier)
E_true = true_efficiency(:,1:T);   % 4 x T (clipped true; for window shaping)
E_ekf  = x_err_ekf(:,1:T);         % 4 x T (signed errors)
E_irls = x_err_irls(:,1:T);        % 4 x T (signed errors)

% 1) Base fault windows from your data-gen schedule
iter_total = size(data.time_arr, 2);
div12 = round(iter_total / 12);
div24 = round(iter_total / 24);
win_mult = [1, 4, 7, 10];   % start multipliers for motors 1..4

% 2) Expand each window to ±3% steady-band (value domain)
%    Window = leave old steady (±3%) -> enter new steady (±3%)
rel_band = 0.03;                 % ±3%
is_fault_window = false(4, T);   % per-motor window mask

for m = 1:4
    % Base window (absolute indices)
    s_abs = win_mult(m) * div12;
    e_abs = s_abs + div24;
    
    % Clamp/map to analysis range [1, T]
    s = max(1, min(T, s_abs));
    e = max(1, min(T, e_abs));
    if s > e, s = e; end
    
    % Neighbor samples to estimate steady levels
    idx_before = max(s-1, 1);
    idx_after  = min(e+1, T);
    E_before = E_true(m, idx_before);
    E_after  = E_true(m, idx_after);
    
    lo_b = (1 - rel_band) * abs(E_before);
    hi_b = (1 + rel_band) * abs(E_before);
    lo_a = (1 - rel_band) * abs(E_after);
    hi_a = (1 + rel_band) * abs(E_after);
    
    % Expand left: back until within old steady band
    s_exp = s;
    while s_exp > 1
        v = abs(E_true(m, s_exp));
        if (v >= lo_b) && (v <= hi_b), break; end
        s_exp = s_exp - 1;
    end
    
    % Expand right: forward until within new steady band
    e_exp = e;
    while e_exp < T
        v = abs(E_true(m, e_exp));
        if (v >= lo_a) && (v <= hi_a), break; end
        e_exp = e_exp + 1;
    end
    
    s_exp = max(1, s_exp);
    e_exp = min(T, e_exp);
    if s_exp > e_exp, s_exp = s; e_exp = e; end
    
    is_fault_window(m, s_exp:e_exp) = true;
end

% 3) Time padding to catch delayed spikes (e.g., EKF lag)
%    (dilate each motor's window by ±PAD_SEC)
PAD_SEC = 0.03; % Adjust if needed
pad = max(1, round(PAD_SEC / dt));
for m = 1:4
    % 1D dilation via convolution
    k = ones(1, 2*pad + 1, 'double');
    dil = conv(double(is_fault_window(m,:)), k, 'same') > 0;
    is_fault_window(m,:) = dil;
end

% 4) Others-only window mask (exclude self)
%    other_event_mask(m,t) = true if any motor ≠ m is in window at t
other_event_mask = false(4, T);
for m = 1:4
    others = setdiff(1:4, m);
    other_event_mask(m, :) = any(is_fault_window(others, :), 1);
end

% Sanity: ensure no self-window contamination
for m = 1:4
    leak = nnz(is_fault_window(m,:) & other_event_mask(m,:));
    if leak > 0
        fprintf('[WARN] Motor %d others-only mask overlaps its own window at %d indices\n', m, leak);
    end
end

% 5) Metrics in others-only windows
%    Spike = max(|x_err_*|); also print the time of that max
rmse_win_ekf  = nan(4,1);
rmse_win_irls = nan(4,1);
std_win_ekf   = nan(4,1);
std_win_irls  = nan(4,1);
max_win_ekf   = nan(4,1);
max_win_irls  = nan(4,1);

for m = 1:4
    mask = other_event_mask(m,:);
    idxs = find(mask);
    if ~isempty(idxs)
        eE = E_ekf(m, idxs);
        eN = E_irls(m, idxs);
        
        % Windowed RMSE/STD (optional)
        rmse_win_ekf(m)  = sqrt(mean(eE.^2));
        rmse_win_irls(m) = sqrt(mean(eN.^2));
        std_win_ekf(m)   = std(eE);
        std_win_irls(m)  = std(eN);
        
        % Spikes
        [max_win_ekf(m), kE]  = max(abs(eE));
        [max_win_irls(m), kN] = max(abs(eN));
        
        % Print when the max spike happened (seconds)
        tE = idxs(kE) * dt;
        tN = idxs(kN) * dt;
        fprintf('Motor %d: others-only samples=%d | EKF max|err|=%.4f @ %.3fs | IRLS max|err|=%.4f @ %.3fs\n', ...
            m, numel(idxs), max_win_ekf(m), tE, max_win_irls(m), tN);
    else
        fprintf('Motor %d: others-only samples=0\n', m);
    end
end

% 6) Plots
metrics = {'RMSE','Std. Dev.','Max Spike'};
vals_irls = [nanmean(rmse_irls), nanmean(std_irls), nanmean(max_win_irls)];
vals_ekf  = [nanmean(rmse_ekf),  nanmean(std_ekf),  nanmean(max_win_ekf)];

fig_b1 = figure('Name','Bars (others-only windows)');
x = 1:numel(metrics);
w = 0.35;
bar(x-w/2, vals_irls, 0.35, 'DisplayName','IRLS');
hold on;
bar(x+w/2, vals_ekf,  0.35, 'DisplayName','EKF');
set(gca,'XTick',x,'XTickLabel',metrics);
ylabel('Value');
grid on;
title('Statistics');
legend('Location','northwest');
exportgraphics(fig_b1, 'bars_event_overall.png', 'Resolution', 300);

fig_b2 = figure('Name','Per-motor max spike');
motors = 1:4;
bar(motors-w/2, max_win_irls, 0.35, 'DisplayName','IRLS'); hold on;
bar(motors+w/2, max_win_ekf,  0.35, 'DisplayName','EKF');
xticks(motors);
xlabel('Motor index');
ylabel('|Max Spike|');
grid on;
title('Per-motor max spike');
legend('Location','northwest');
exportgraphics(fig_b2, 'bars_event_permotor_max.png', 'Resolution', 300);

%==================%
% Print statistics %
%==================%
fprintf("\n=== EKF stats ===\n");
fprintf("RMSE       = (%f, %f, %f, %f)\n", rmse_ekf(1), rmse_ekf(2), rmse_ekf(3), rmse_ekf(4));
fprintf("Variance   = (%f, %f, %f, %f)\n", var_ekf(1),  var_ekf(2),  var_ekf(3),  var_ekf(4));
fprintf("Std. Dev   = (%f, %f, %f, %f)\n", std_ekf(1),  std_ekf(2),  std_ekf(3),  std_ekf(4));
fprintf("Max spike  = (%f, %f, %f, %f)\n", max_win_ekf(1), max_win_ekf(2), max_win_ekf(3), max_win_ekf(4))

fprintf("\n=== IRLS stats ===\n");
fprintf("RMSE       = (%f, %f, %f, %f)\n", rmse_irls(1), rmse_irls(2), rmse_irls(3), rmse_irls(4));
fprintf("Variance   = (%f, %f, %f, %f)\n", var_irls(1),  var_irls(2),  var_irls(3),  var_irls(4));
fprintf("Std. Dev   = (%f, %f, %f, %f)\n", std_irls(1),  std_irls(2),  std_irls(3),  std_irls(4));
fprintf("Max spike  = (%f, %f, %f, %f)\n", max_win_irls(1), max_win_irls(2), max_win_irls(3), max_win_irls(4))

fprintf('\n=== Averaging ===\n');
fprintf('RMSE : IRLS=%.4f, EKF=%.4f\n', nanmean(rmse_irls), nanmean(rmse_ekf));
fprintf('Max  : IRLS=%.4f, EKF=%.4f\n\n',  nanmean(max_win_irls),  nanmean(max_win_ekf));

disp("Press any key to leave");
pause;
close all;

%=======================%
% Helper: get_new_batch %
%=======================%
function batch = get_new_batch(data, iteration, len)
idx = iteration;
idx_end = idx + len - 1;
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