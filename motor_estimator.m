classdef motor_estimator
    properties
        math;
        
        n = 10; % Number of trajectory points (-1 for # of segments)
        I = eye(3);
        g = 9.8;
        e3 = [0; 0; 1];
        dt;
        dt2; % dt squared
        
        G_i = ...
            [60, 0, 0, 0, 0, 0, 0, 0, 0, 0;  % G_vx
            0, 60, 0, 0, 0, 0, 0, 0, 0, 0;  % G_vy
            0, 0, 60, 0, 0, 0, 0, 0, 0, 0;  % G_vz
            0, 0, 0, 30, 0, 0, 0, 0, 0, 0;  % G_px
            0, 0, 0, 0, 30, 0, 0, 0, 0, 0;  % G_py
            0, 0, 0, 0, 0, 30, 0, 0, 0, 0;  % G_pz
            0, 0, 0, 0, 0, 0, 1, 0, 0, 0;  % G_Wx
            0, 0, 0, 0, 0, 0, 0, 1, 0, 0;  % G_Wy
            0, 0, 0, 0, 0, 0, 0, 0, 1, 0;  % G_Wz
            0, 0, 0, 0, 0, 0, 0, 0, 0, 1]; % G_R
        G = [];
        
        Dphi = [1, 0, 0, 0; ...
            0, 1, 0, 0; ...
            0, 0, 1, 0; ...
            0, 0, 0, 1; ...
            -1, 0, 0, 0; ...
            0, -1, 0, 0; ...
            0, 0, -1, 0; ...
            0, 0, 0, -1];
        Dphi_t;
        
        mu_smooth = 5e-3;
        
        lower_bound = 0;
        upper_bound = 1;
        bound_safe_eps = 1e-3;
        
        cond_threshold = 1e4;
        rcond_threshold;
    end
    
    methods
        function ret = init(obj, math, dt)
            obj.math = math;
            obj.dt = dt;
            obj.dt2 = dt^2;
            obj.rcond_threshold = 1 / obj.cond_threshold;
            
            for i = 1:obj.n-1
                obj.G = blkdiag(obj.G, obj.G_i);
            end
            
            obj.Dphi_t = obj.Dphi.';
            
            ret = obj;
        end
        
        function ret = calc_trajectory_residual_vector(obj, m, J, x, f_motors, v, p, W, R)
            J_inv = diag([1/J(1,1); 1/J(2,2); 1/J(3,3)]);
            residual = zeros(10 * (obj.n - 1), 1);
            
            % FIXME
            d = 0.225;
            c = 0.009012;
            
            for i=1:obj.n-1
                f1 = f_motors(1, i);
                f2 = f_motors(2, i);
                f3 = f_motors(3, i);
                f4 = f_motors(4, i);
                
                f = x(1)*f1 + ...
                    x(2)*f2 + ...
                    x(3)*f3 + ...
                    x(4)*f4;
                M = [d*(-x(1)*f1 + x(2)*f2 + x(3)*f3 - x(4)*f4); ...
                    d*(+x(1)*f1 + x(2)*f2 - x(3)*f3 - x(4)*f4); ...
                    c*(-x(1)*f1 + x(2)*f2 - x(3)*f3 + x(4)*f4)];
                
                ge3_fRe3m = obj.g*obj.e3 - (f*R(:, :, i)*obj.e3)/m;
                pred_v = v(:, i) + ge3_fRe3m*obj.dt;
                pred_p = p(:, i) + v(:, i)*obj.dt + ge3_fRe3m*obj.dt2/2;
                
                Jinv_M_WJW = J_inv * (M - cross(W(:, i), J*W(:, i)));
                pred_W = W(:, i) + Jinv_M_WJW*obj.dt;
                pred_R = R(:, :, i) * (obj.I + obj.math.hat_map_3x3(W(:, i))*obj.dt);
                %pred_R = R(:, :, i) * expm(obj.math.hat_map_3x3(W(:, i)) * obj.dt);
                
                delta_R = R(:, :, i).' * R(:, :, i+1);
                delta_pred_R = pred_R.' * R(:, :, i+1);
                
                % Velocity error (1~3)
                residual(10*(i-1)+1:10*(i-1)+3) = v(:, i+1) - pred_v;
                % Position error (4~6)
                residual(10*(i-1)+4:10*(i-1)+6) = p(:, i+1) - pred_p;
                % Angular velocity error (7~9)
                residual(10*(i-1)+7:10*(i-1)+9) = W(:, i+1) - pred_W;
                % Rotational error (10)
                residual(10*(i-1)+10:10*(i-1)+10) = trace(obj.I - delta_R.' * delta_pred_R) / 2;
            end
            
            ret = residual / sqrt(obj.n - 1);
        end
        
        function phi = calc_phi_vector(obj, x)
            lb = obj.lower_bound - obj.bound_safe_eps;
            hb = obj.upper_bound + obj.bound_safe_eps;
            phi = [x(1) - hb; ...
                x(2) - hb; ...
                x(3) - hb; ...
                x(4) - hb; ...
                lb - x(1); ...
                lb - x(2); ...
                lb - x(3); ...
                lb - x(4)];
        end
        
        function [r_dual, r_cent] = calc_primal_dual_residual_vector(obj, x, x_last, phi, Jt, f_residual, lambda, t)
            grad_smooth = obj.mu_smooth * (x - x_last);
            r_dual = Jt*obj.G*f_residual + obj.Dphi_t*lambda + grad_smooth;
            r_cent = -diag(lambda)*phi - (1/t)*ones(8, 1);
        end
        
        function KKT = calc_primal_dual_kkt_matrix(obj, x, x_last, phi, Jt_G_J, lambda)
            H_smooth = obj.mu_smooth * eye(4);
            KKT = [Jt_G_J + H_smooth, obj.Dphi_t; ...
                -diag(lambda)*obj.Dphi, -diag(phi)];
        end
        
        function gap = calc_duality_gap(obj, phi, lambda)
            gap = -phi.' * lambda;
        end
        
        function ret = calc_trajectory_jacobian(obj, m, J, c, d, f_motors, v, p, W, R)
            X_DIM = 4; % [eta1; eta2; eta3; et4]
            RES_DIM = 10; % [res_v, res_p, res_W, res_R]
            Jf = zeros(RES_DIM*(obj.n - 1), X_DIM);
            
            Jx = J(1,1);
            Jy = J(2,2);
            Jz = J(3,3);
            
            for i=1:obj.n-1
                f1 = f_motors(1, i);
                f2 = f_motors(2, i);
                f3 = f_motors(3, i);
                f4 = f_motors(4, i);
                
                delta_R = R(:,:,i).'*R(:,:,i+1);
                r31 = R(3,1,i);
                r32 = R(3,2,i);
                r33 = R(3,3,i);
                
                fR = [f1*r31, f2*r31, f3*r31, f4*r31; ...
                    f1*r32, f2*r32, f3*r32, f4*r32; ...
                    f1*r33, f2*r33, f3*r33, f4*r33];
                
                % d(res_v)/dx
                r_start = 10*(i-1)+1;
                r_end =   10*(i-1)+3;
                Jf(r_start:r_end, 1:4) = obj.dt/m * fR;
                
                % d(res_p)/dx
                r_start = 10*(i-1)+4;
                r_end =   10*(i-1)+6;
                Jf(r_start:r_end, 1:4) = obj.dt2/(2*m) * fR;
                
                % d(res_W)/dx
                r_start = 10*(i-1)+7;
                r_end =   10*(i-1)+9;
                Jf(r_start:r_end,1:4) = obj.dt * ...
                    [+f1*d/Jx, -f2*d/Jx, -f3*d/Jx, +f4*d/Jx; ...
                    -f1*d/Jy, -f2*d/Jy, +f3*d/Jy, +f4*d/Jy; ...
                    +f1*c/Jz, -f2*c/Jz, +f3*c/Jz, -f4*c/Jz];
                
                % d(res_R)/dx
                r_start = 10*(i-1)+10;
                r_end =   10*(i-1)+10;
                r23_r32 = delta_R(2,3)-delta_R(3,2);
                r31_r13 = delta_R(3,1)-delta_R(1,3);
                r12_r21 = delta_R(1,2)-delta_R(2,1);
                r1 = -r23_r32*f1*d/Jx ...
                    +r31_r13*f1*d/Jy ...
                    -r12_r21*f1*c/Jz;
                r2 = +r23_r32*f2*d/Jx ...
                    +r31_r13*f2*d/Jy ...
                    +r12_r21*f2*c/Jz;
                r3 = +r23_r32*f3*d/Jx ...
                    -r31_r13*f3*d/Jy ...
                    -r12_r21*f3*c/Jz;
                r4 = -r23_r32*f4*d/Jx ...
                    -r31_r13*f4*d/Jy ...
                    +r12_r21*f4*c/Jz;
                Jf(r_start:r_end, 1:4) = obj.dt2 * [r1, r2, r3, r4] / 2;
            end
            
            ret = Jf;
        end
        
        function [ret_x, cond] = run_primal_dual(obj, iteration, batch, x, x_last)
            x_arr = [];
            iteration_arr = [];
            residual_arr = [];
            r_dual_now_arr = [];
            gap_arr = [];
            iteration = 0;
            
            iteration = iteration + 1;
            iteration_arr(end+1) = iteration;
            x_arr(:, end+1) = x;
            
            lb = obj.lower_bound - obj.bound_safe_eps;
            hb = obj.upper_bound + obj.bound_safe_eps;
            y = zeros(12, 1); % x (4) + lambda (8)
            y(1:4) = x;
            y(5:12) = [1 / (hb - x(1)); ...
                1 / (hb - x(2)); ...
                1 / (hb - x(3)); ...
                1 / (hb - x(4)); ...
                1 / (x(1) - lb); ...
                1 / (x(2) - lb); ...
                1 / (x(3) - lb); ...
                1 / (x(4) - lb)];
            lambda = y(5:12);
            
            % Linearization
            Jf = calc_trajectory_jacobian(obj, batch.m, batch.J, batch.c, batch.d, ...
                batch.f_motors, batch.v, batch.p, batch.W, batch.R);
            Jf_t = Jf.';
            Jt_G_J = Jf_t * obj.G * Jf;
            
            % Skip Degenerate Jacobian due to insufficient excitement of the trajectory
            rcond_JtJ = rcond(Jf.'*obj.G*Jf);
            cond = 1 / rcond_JtJ;
            if rcond_JtJ < obj.rcond_threshold
                fprintf('x: Degenerate Jacobian detected – skip this time step\n');
                ret_x = x_last;
                return;
            end
            
            alpha = 0.1; % typically 0.01 to 0.1
            beta = 0.5; % typically 0.3 to 0.8
            mu = 4.5;
            m = 8; % Constraints numbers (i.e., lower_bound < x < upper_bound)
            
            for i = 1:200 % Limit on maximum iteration
                %fprintf("iteration: %d\n", iteration)
                
                % Construct constraint vector
                phi = calc_phi_vector(obj, x);
                
                % Deterimine t size of the log barrier
                gap = calc_duality_gap(obj, phi, lambda);
                t = mu * m / gap;
                
                % Solve primal dual step
                f_residual = calc_trajectory_residual_vector(obj, batch.m, batch.J, x, batch.f_motors, batch.v, batch.p, batch.W, batch.R);
                [r_dual, r_cent] = calc_primal_dual_residual_vector(obj, x, x_last, phi, Jf_t, f_residual, lambda, t);
                pd_residual = [r_dual; r_cent];
                KKT = calc_primal_dual_kkt_matrix(obj, x, x_last, phi, Jt_G_J, lambda);
                delta_y = -KKT \ pd_residual;
                
                % Calculate maximum step size that won't violate KKT conditions
                delta_lambda = delta_y(5:12);
                idx = delta_lambda < 0;
                if any(idx)
                    ratios = -lambda(idx) ./ delta_lambda(idx);
                    s_max = min(1, min(ratios));
                else
                    s_max = 1;
                end
                s = 0.99 * s_max; % Prevent singluar value at the edge of the constraints
                
                if 0
                    % Full step
                    y = y + s*delta_y;
                    phi = calc_phi_vector(obj, y(1:4));
                    f_residual = calc_trajectory_residual_vector(obj, batch.m, batch.J, y(1:4), batch.f_motors, batch.v, batch.p, batch.W, batch.R);
                    [r_dual_now, r_cent_now] = calc_primal_dual_residual_vector(obj, x, x_last, phi, Jf_t, f_residual, y(5:12), t);
                    r_now = [r_dual_now; r_cent_now];
                else
                    % Line backtracking to prevent dual variable to be negative
                    f_residual = calc_trajectory_residual_vector(obj, batch.m, batch.J, x, batch.f_motors, batch.v, batch.p, batch.W, batch.R);
                    [r_dual_last, r_cent_last] = calc_primal_dual_residual_vector(obj, x, x_last, phi, Jf_t, f_residual, lambda, t);
                    r_last = [r_dual_last; r_cent_last];
                    y0 = y;
                    iter = 0;
                    while iter < 50
                        % Calculate backtracking residual
                        y = y0 + s*delta_y;
                        phi = calc_phi_vector(obj, y(1:4));
                        f_residual = calc_trajectory_residual_vector(obj, batch.m, batch.J, y(1:4), batch.f_motors, batch.v, batch.p, batch.W, batch.R);
                        [r_dual_now, r_cent_now] = calc_primal_dual_residual_vector(obj, x, x_last, phi, Jf_t, f_residual, y(5:12), t);
                        r_now = [r_dual_now; r_cent_now];
                        
                        % Stop if Armijo condition is fufilled
                        if norm(r_now) < (1 - alpha*s)*norm(r_last) + 1e-3 % Add tolerance range due to noise
                            break;
                        end
                        
                        % Update step size
                        s = s * beta;
                        
                        % Update iteration
                        iter = iter + 1;
                    end
                end
                
                % Update primal and dual variables
                x = y(1:4);
                lambda = y(5:12);
                gap = calc_duality_gap(obj, phi, lambda);
                
                % Profiling
                iteration = iteration + 1;
                iteration_arr(end+1) = iteration;
                x_arr(:, end+1) = x;
                residual_arr(end+1) = sqrt(f_residual.' * obj.G * f_residual);
                r_dual_now_arr(end+1) = norm(r_dual_now);
                gap_arr(end+1) = gap;
                
                if norm(r_dual_now) < 1e-6 && gap < 1e-6
                    break;
                end
            end
            
            % Profiling
            if 1
                fig = figure('Name', 'Efficiency vs Iteration');
                subplot (4, 1, 1);
                h1 = plot(iteration_arr - 1, x_arr(1, :), 'LineWidth', 1.7);
                h2 = yline(batch.motor_efficiency(1), '--k', 'Color', 'r', 'LineWidth', 1.7);
                legend('Estimated', 'True', 'Location', 'southeast', 'Orientation', 'horizontal');
                title('Motor efficiency: Estimated vs True', 'FontSize', 11);
                xlabel('Iteration number', 'FontSize', 11);
                ylabel('\eta_1', 'FontSize', 11);
                xlim([0, iteration - 1]);
                ylim([0.2 1.2]);
                grid on;
                subplot (4, 1, 2);
                h1 = plot(iteration_arr - 1, x_arr(2, :), 'LineWidth', 1.7);
                h2 = yline(batch.motor_efficiency(2), '--k', 'Color', 'r', 'LineWidth', 1.7);
                legend('Estimated', 'True', 'Location', 'southeast', 'Orientation', 'horizontal');
                xlabel('Iterations', 'FontSize', 11);
                ylabel('\eta_2', 'FontSize', 11);
                xlim([0, iteration - 1]);
                ylim([0.2 1.2]);
                grid on;
                subplot (4, 1, 3);
                h1 = plot(iteration_arr - 1, x_arr(3, :), 'LineWidth', 1.7);
                h2 = yline(batch.motor_efficiency(3), '--k', 'Color', 'r', 'LineWidth', 1.7);
                legend('Estimated', 'True', 'Location', 'southeast', 'Orientation', 'horizontal');
                xlabel('Iterations', 'FontSize', 11);
                ylabel('\eta_3', 'FontSize', 11);
                xlim([0, iteration - 1]);
                ylim([0.2 1.2]);
                grid on;
                subplot (4, 1, 4);
                h1 = plot(iteration_arr - 1, x_arr(4, :), 'LineWidth', 1.7);
                h2 = yline(batch.motor_efficiency(4), '--k', 'Color', 'r', 'LineWidth', 1.7);
                legend('Estimated', 'True', 'Location', 'southeast', 'Orientation', 'horizontal');
                xlabel('Iterations', 'FontSize', 11);
                ylabel('\eta_4', 'FontSize', 11);
                xlim([0, iteration - 1]);
                ylim([0.2 1.2]);
                grid on;
                exportgraphics(fig, 'efficiency_iteration.png');
                
                fig = figure('Name', 'Residual vs Iteration');
                plot(iteration_arr(2:end) - 1, residual_arr(1:end), 'o-', ...
                    'LineWidth', 1.7, ...
                    'Color', '#1f77b4', ...
                    'MarkerEdgeColor', '#1f77b4', ...
                    'MarkerFaceColor', '#1f77b4');
                xlabel('Iterations', 'FontSize', 11);
                ylabel('Residual norm', 'FontSize', 11);
                xlim([1, iteration - 1]);
                title('Residual norm vs Iteration', 'FontSize', 11);
                grid on;
                exportgraphics(fig, 'residual_iteration.png');
                
                fig = figure('Name', 'KKT dual residual vs Iteration');
                plot(iteration_arr(2:end) - 1, r_dual_now_arr(1:end), 'o-', ...
                    'LineWidth', 1.7, ...
                    'Color', '#1f77b4', ...
                    'MarkerEdgeColor', '#1f77b4', ...
                    'MarkerFaceColor', '#1f77b4');
                xlabel('Iterations', 'FontSize', 11);
                ylabel('$\|\textbf{r}_{\textbf{dual}}\|$', 'Interpreter', 'latex', 'FontSize', 11);
                xlim([1, iteration - 1]);
                title('KKT dual residual convergence', 'FontSize', 11);
                grid on;
                
                fig = figure('Name', 'Gap vs Iteration');
                plot(iteration_arr(2:end) - 1, gap_arr(1:end), 'o-', ...
                    'LineWidth', 1.7, ...
                    'Color', '#1f77b4', ...
                    'MarkerEdgeColor', '#1f77b4', ...
                    'MarkerFaceColor', '#1f77b4');
                xlabel('Iterations', 'FontSize', 11);
                ylabel('$\mathbf{\hat{\delta}}$', 'Interpreter', 'latex', 'FontSize', 11);
                xlim([1, iteration - 1]);
                title('Surrogate duality gap convergence', 'FontSize', 11);
                grid on;
                
                fig = figure('Name', 'KKT dual residual norm vs Iteration');
                subplot(2, 1, 1);
                plot(iteration_arr(2:end) - 1, r_dual_now_arr(1:end), 'o-', ...
                    'LineWidth', 1.7, ...
                    'Color', '#1f77b4', ...
                    'MarkerEdgeColor', '#1f77b4', ...
                    'MarkerFaceColor', '#1f77b4');
                xlabel('Iterations', 'FontSize', 11);
                ylabel('$\|\textbf{r}_{\textbf{dual}}\|$', 'Interpreter', 'latex', 'FontSize', 11);
                xlim([1, iteration - 1]);
                title('KKT dual residual convergence (log)', 'FontSize', 11);
                grid on;
                set(gca, 'YScale', 'log');
                yticks(10.^(-10:1));
                %
                subplot(2, 1, 2);
                plot(iteration_arr(2:end) - 1, gap_arr(1:end), 'o-', ...
                    'LineWidth', 1.7, ...
                    'Color', '#1f77b4', ...
                    'MarkerEdgeColor', '#1f77b4', ...
                    'MarkerFaceColor', '#1f77b4');
                xlabel('Iterations', 'FontSize', 11);
                ylabel('$\mathbf{\hat{\delta}}$', 'Interpreter', 'latex', 'FontSize', 11);
                xlim([1, iteration - 1]);
                title('Surrogate duality gap convergence  (log)', 'FontSize', 11);
                grid on;
                set(gca, 'YScale', 'log');
                yticks([1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1]);
                exportgraphics(fig, 'kkt_convergence.png');
                
                pause;
                close all;
            end
            
            ret_x = x;
        end
        
        function [ret_x, cond_max] = run(obj, iteration, batch, x, x_last)
            [ret_x, cond_max] = run_primal_dual(obj, iteration, batch, x, x_last);
        end
    end
end