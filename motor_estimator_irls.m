classdef motor_estimator_irls
    properties
        math;
        
        ENABLE_PLOTS = false;
        
        n  = 10;
        I  = eye(3);
        g  = 9.8;
        e3 = [0; 0; 1];
        dt;
        dt2; % dt^2
        
        G_i = [50, 0, 0, 0, 0, 0, 0, 0, 0, 0;  % v_x
            0, 50, 0, 0, 0, 0, 0, 0, 0, 0;     % v_y
            0, 0, 50, 0, 0, 0, 0, 0, 0, 0;     % v_z
            0, 0, 0, 50, 0, 0, 0, 0, 0, 0;     % p_x
            0, 0, 0, 0, 50, 0, 0, 0, 0, 0;     % p_y
            0, 0, 0, 0, 0, 50, 0, 0, 0, 0;     % p_z
            0, 0, 0, 0, 0, 0, 1, 0, 0, 0;      % W_x
            0, 0, 0, 0, 0, 0, 0, 1, 0, 0;      % W_y
            0, 0, 0, 0, 0, 0, 0, 0, 1, 0;      % W_z
            0, 0, 0, 0, 0, 0, 0, 0, 0, 50];    % R
        G = [];
        
        % Constraints：0 <= eta <= 1
        Dphi = [1, 0, 0, 0;
            0, 1, 0, 0;
            0, 0, 1, 0;
            0, 0, 0, 1;
            -1, 0, 0, 0;
            0,-1, 0, 0;
            0, 0,-1, 0;
            0, 0, 0,-1];
        Dphi_t;
        
        gamma_smooth = 5e-3;
        
        lower_bound = 0;
        upper_bound = 1;
        bound_safe_eps = 1e-3;
        
        cond_threshold  = 1e4;
        rcond_threshold;
        
        % Robust z-score with median absolute deviation (MAD)
        use_mad_filter = true;
        z_soft         = 2.5;    % Soft drop enabling
        mad_soft_power = 7.0;    % Soft-rejection power
        w_min          = 0.03;   % Min weight
        mad_hard_drop  = true;   % Hard drop enabling
        z_hard         = 2.5;    % Hard-rejection
        mad_eps        = 1e-9;   % Prevent div by zero
    end
    
    methods
        function ret = init(obj, math, dt)
            obj.math = math;
            obj.dt   = dt;
            obj.dt2  = dt^2;
            obj.rcond_threshold = 1 / obj.cond_threshold;
            
            % Initial G
            obj.G = [];
            for i = 1:obj.n-1
                obj.G = blkdiag(obj.G, obj.G_i);
            end
            
            obj.Dphi_t = obj.Dphi.';
            ret = obj;
        end
        
        function ret = calc_trajectory_residual_vector(obj, m, J, c, d, x, f_motors, v, p, W, R)
            J_inv = diag([1/J(1,1); 1/J(2,2); 1/J(3,3)]);
            residual = zeros(10 * (obj.n - 1), 1);
            
            for i = 1:obj.n-1
                f1 = f_motors(1, i);
                f2 = f_motors(2, i);
                f3 = f_motors(3, i);
                f4 = f_motors(4, i);
                
                f = x(1)*f1 + x(2)*f2 + x(3)*f3 + x(4)*f4;
                M = [ d*(-x(1)*f1 + x(2)*f2 + x(3)*f3 - x(4)*f4);
                    d*(+x(1)*f1 + x(2)*f2 - x(3)*f3 - x(4)*f4);
                    c*(-x(1)*f1 + x(2)*f2 - x(3)*f3 + x(4)*f4) ];
                
                ge3_fRe3m = obj.g*obj.e3 - (f*R(:,:,i)*obj.e3)/m;
                pred_v = v(:,i) + ge3_fRe3m*obj.dt;
                pred_p = p(:,i) + v(:,i)*obj.dt + ge3_fRe3m*(obj.dt2/2);
                
                Jinv_M_WJW = J_inv * (M - cross(W(:,i), J*W(:,i)));
                pred_W = W(:,i) + Jinv_M_WJW*obj.dt;
                pred_R = R(:,:,i) * (obj.I + obj.math.hat_map_3x3(pred_W * obj.dt));
                % pred_R = R(:,:,i) * expm(obj.math.hat_map_3x3(pred_W * obj.dt));
                
                delta_R = R(:,:,i).' * R(:,:,i+1);
                delta_pred_R = R(:,:,i).' * pred_R;
                
                residual(10*(i-1)+1:10*(i-1)+3)  = v(:,i+1) - pred_v; % v
                residual(10*(i-1)+4:10*(i-1)+6)  = p(:,i+1) - pred_p; % p
                residual(10*(i-1)+7:10*(i-1)+9)  = W(:,i+1) - pred_W; % W
                residual(10*(i-1)+10:10*(i-1)+10)= trace(obj.I - delta_R.'*delta_pred_R)/2; % R
            end
            ret = residual;
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
            grad_smooth = obj.gamma_smooth * (x - x_last);
            r_dual = Jt * obj.G * f_residual + obj.Dphi_t * lambda + grad_smooth;
            r_cent = -diag(lambda) * phi - (1/t) * ones(8,1);
        end
        
        function KKT = calc_primal_dual_kkt_matrix(obj, phi, Jt_G_J, lambda)
            H_smooth = obj.gamma_smooth * eye(4);
            KKT = [Jt_G_J + H_smooth, obj.Dphi_t;
                -diag(lambda)*obj.Dphi, -diag(phi)];
        end
        
        function gap = calc_duality_gap(obj, phi, lambda)
            gap = -phi.' * lambda;
        end
        
        function Jf = calc_trajectory_jacobian(obj, m, J, c, d, f_motors, R)
            X_DIM = 4; RES_DIM = 10;
            Jf = zeros(RES_DIM*(obj.n - 1), X_DIM);
            
            Jx = J(1,1); Jy = J(2,2); Jz = J(3,3);
            
            for i=1:obj.n-1
                f1 = f_motors(1,i);
                f2 = f_motors(2,i);
                f3 = f_motors(3,i);
                f4 = f_motors(4,i);
                
                delta_R = R(:,:,i).' * R(:,:,i+1);
                r13 = R(1,3,i); r23 = R(2,3,i); r33 = R(3,3,i);
                
                fR = [f1*r13, f2*r13, f3*r13, f4*r13;
                    f1*r23, f2*r23, f3*r23, f4*r23;
                    f1*r33, f2*r33, f3*r33, f4*r33];
                
                % d(res_v)/dx
                r_start = 10*(i-1)+1; r_end = 10*(i-1)+3;
                Jf(r_start:r_end,1:4) = obj.dt/m * fR;
                
                % d(res_p)/dx
                r_start = 10*(i-1)+4; r_end = 10*(i-1)+6;
                Jf(r_start:r_end,1:4) = (obj.dt2/(2*m)) * fR;
                
                % d(res_W)/dx
                r_start = 10*(i-1)+7; r_end = 10*(i-1)+9;
                Jf(r_start:r_end,1:4) = obj.dt * [ +f1*d/Jx, -f2*d/Jx, -f3*d/Jx, +f4*d/Jx;
                    -f1*d/Jy, -f2*d/Jy, +f3*d/Jy, +f4*d/Jy;
                    +f1*c/Jz, -f2*c/Jz, +f3*c/Jz, -f4*c/Jz ];
                
                % d(res_R)/dx
                r_start = 10*(i-1)+10; r_end = 10*(i-1)+10;
                r23_r32 = delta_R(2,3)-delta_R(3,2);
                r31_r13 = delta_R(3,1)-delta_R(1,3);
                r12_r21 = delta_R(1,2)-delta_R(2,1);
                r1 = -r23_r32*f1*d/Jx + r31_r13*f1*d/Jy - r12_r21*f1*c/Jz;
                r2 = +r23_r32*f2*d/Jx + r31_r13*f2*d/Jy + r12_r21*f2*c/Jz;
                r3 = +r23_r32*f3*d/Jx - r31_r13*f3*d/Jy - r12_r21*f3*c/Jz;
                r4 = -r23_r32*f4*d/Jx - r31_r13*f4*d/Jy + r12_r21*f4*c/Jz;
                Jf(r_start:r_end,1:4) = (obj.dt2/2) * [r1, r2, r3, r4];
            end
        end
        
        function [w_seg, z, e_seg, med_e, mad_e] = mad_outlier_weights(obj, f_residual)
            % Calculate weighted residual
            m1 = obj.n - 1;
            e_seg = zeros(m1,1);
            for i = 1:m1
                idx = 10*(i-1)+1 : 10*i;
                ri  = f_residual(idx);
                e_seg(i) = ri.' * obj.G_i * ri;
            end
            
            % Calculate robust z-score with median absolute deviation (MAD)
            med_e     = median(e_seg);
            mad_e_raw = median(abs(e_seg - med_e));
            mad_e     = 1.4826 * mad_e_raw;
            
            if mad_e < obj.mad_eps
                w_seg = ones(m1,1); z = zeros(m1,1);
                return;
            end
            
            % Calculate new weights
            z = abs(e_seg - med_e) ./ max(mad_e, obj.mad_eps);
            w_seg = 1 ./ (1 + (z./obj.z_soft).^obj.mad_soft_power);
            
            if obj.mad_hard_drop
                w_seg(z > obj.z_hard) = 0;
            else
                w_seg = max(w_seg, obj.w_min);
            end
            
            medw = median(w_seg);
            if medw > 0, w_seg = w_seg / medw; end
        end
        
        function [ret_x, cond_out] = run_primal_dual(obj, batch, x, x_last)
            iteration      = 0;
            iteration_arr  = [];
            x_arr          = [];
            residual_arr   = [];
            residual_per_k = {};
            r_dual_now_arr = [];
            gap_arr        = [];
            
            % Compute Jacobian matrix (independent from states, compute
            % once)
            Jf   = obj.calc_trajectory_jacobian(batch.m, batch.J, batch.c, batch.d, ...
                batch.f_motors, batch.R);
            Jf_t = Jf.';
            
            % Line search and barrier method's parameters
            alpha = 0.1;
            beta = 0.5;
            mu = 4.0;
            m_bar = 8;
            max_inner = 80;
            max_irls  = 3;
            
            x_cur = x;
            cond_out = NaN;
            
            % Outer loop: IRLS
            for ir = 1:max_irls
                residual_this_round = [];
                
                f_res0 = obj.calc_trajectory_residual_vector( ...
                    batch.m, batch.J, batch.c, batch.d, ...
                    x_cur, batch.f_motors, batch.v, batch.p, batch.W, batch.R);
                
                if obj.use_mad_filter
                    [w_seg, ~, ~, ~, ~] = obj.mad_outlier_weights(f_res0);
                else
                    w_seg = ones(obj.n-1,1);
                end
                
                Gblk = [];
                for ii = 1:(obj.n-1)
                    Gblk = blkdiag(Gblk, w_seg(ii) * obj.G_i);
                end
                obj.G = Gblk;
                
                Jt_G_J  = Jf_t * obj.G * Jf;
                
                rcond_JtJ = rcond(Jt_G_J);
                cond_out  = 1 / rcond_JtJ;
                if rcond_JtJ < obj.rcond_threshold
                    fprintf('x: Degenerate Jacobian (IRLS %d) – skip this time step\n', ir);
                    ret_x = x_last;
                    return;
                end
                
                % Inner loop: interior point method
                lb = obj.lower_bound - obj.bound_safe_eps;
                hb = obj.upper_bound + obj.bound_safe_eps;
                
                y = zeros(12,1);
                y(1:4) = x_cur;
                y(5:12) = [1/(hb-x_cur(1)); ...
                    1/(hb-x_cur(2)); ...
                    1/(hb-x_cur(3)); ...
                    1/(hb-x_cur(4)); ...
                    1/(x_cur(1)-lb); ...
                    1/(x_cur(2)-lb); ...
                    1/(x_cur(3)-lb); ...
                    1/(x_cur(4)-lb)];
                lambda = y(5:12);
                
                for it = 1:max_inner
                    xk  = y(1:4);
                    phi = obj.calc_phi_vector(xk);
                    gap = obj.calc_duality_gap(phi, lambda);
                    t = mu * m_bar / max(gap, 1e-12);
                    
                    f_res = obj.calc_trajectory_residual_vector( ...
                        batch.m, batch.J, batch.c, batch.d, ...
                        xk, batch.f_motors, batch.v, batch.p, batch.W, batch.R);
                    
                    [r_dual, r_cent] = obj.calc_primal_dual_residual_vector( ...
                        xk, x_last, phi, Jf_t, f_res, lambda, t);
                    KKT   = obj.calc_primal_dual_kkt_matrix(phi, Jt_G_J, lambda);
                    rhs   = -[r_dual; r_cent];
                    delta_y = KKT \ rhs;
                    
                    % Max step size
                    delta_lambda = delta_y(5:12);
                    idx = (delta_lambda < 0);
                    if any(idx)
                        s_max = min(1, min(-lambda(idx)./delta_lambda(idx)));
                    else
                        s_max = 1;
                    end
                    s = 0.99 * s_max;
                    
                    % Armijo backtracking
                    y0 = y;
                    r_last = [r_dual; r_cent];
                    r_now  = r_last;
                    while true
                        y_try   = y0 + s * delta_y;
                        x_try   = y_try(1:4);
                        lam_try = y_try(5:12);
                        
                        phi_try = obj.calc_phi_vector(x_try);
                        f_res_t = obj.calc_trajectory_residual_vector( ...
                            batch.m, batch.J, batch.c, batch.d, ...
                            x_try, batch.f_motors, batch.v, batch.p, batch.W, batch.R);
                        [r_dual_now, r_cent_now] = obj.calc_primal_dual_residual_vector( ...
                            x_try, x_last, phi_try, Jf_t, f_res_t, lam_try, t);
                        r_now = [r_dual_now; r_cent_now];
                        
                        if norm(r_now) <= (1 - alpha*s) * norm(r_last) + 1e-3 || s < 1e-6
                            y = y_try;
                            lambda = lam_try;
                            
                            iteration = iteration + 1;
                            iteration_arr(end+1) = iteration;
                            x_arr(:, end+1) = y(1:4);
                            residual_arr(end+1) = sqrt(f_res_t.' * obj.G * f_res_t);
                            residual_this_round(end+1) = sqrt(f_res_t.' * obj.G * f_res_t);
                            r_dual_now_arr(end+1) = norm(r_dual_now);
                            gap_arr(end+1) = obj.calc_duality_gap(phi_try, lam_try);
                            
                            break;
                        end
                        s = beta * s;
                    end
                    
                    % Inner loop stop condition
                    if norm(r_now) < 1e-6 && obj.calc_duality_gap(phi, lambda) < 1e-6
                        break;
                    end
                end
                
                x_new = y(1:4);
                
                % Logging k-th outer-loop residuals
                residual_per_k{end+1} = residual_this_round;
                
                % Outer loop stop condition
                if norm(x_new - x_cur, inf) < 1e-3
                    x_cur = x_new;
                    break;
                end
                
                x_cur = x_new;
            end
            
            % Profiling
            if obj.ENABLE_PLOTS
                fig = figure('Name', 'Efficiency vs Iteration');
                subplot (4, 1, 1);
                h1 = plot(iteration_arr - 1, x_arr(1, :), 'LineWidth', 2);
                h2 = yline(batch.motor_efficiency(1), '--k', 'Color', 'r', 'LineWidth', 2);
                legend('Estimated', 'True', 'Location', 'southeast', 'Orientation', 'horizontal');
                title('Motor efficiency: Estimated vs True', 'FontSize', 11);
                ylabel('\eta_1', 'FontSize', 11);
                xlim([0, iteration - 1]);
                ylim([0.2 1.2]);
                grid on;
                %
                subplot (4, 1, 2);
                h1 = plot(iteration_arr - 1, x_arr(2, :), 'LineWidth', 2);
                h2 = yline(batch.motor_efficiency(2), '--k', 'Color', 'r', 'LineWidth', 2);
                legend('Estimated', 'True', 'Location', 'southeast', 'Orientation', 'horizontal');
                ylabel('\eta_2', 'FontSize', 11);
                xlim([0, iteration - 1]);
                ylim([0.2 1.2]);
                grid on;
                %
                subplot (4, 1, 3);
                h1 = plot(iteration_arr - 1, x_arr(3, :), 'LineWidth', 2);
                h2 = yline(batch.motor_efficiency(3), '--k', 'Color', 'r', 'LineWidth', 2);
                legend('Estimated', 'True', 'Location', 'southeast', 'Orientation', 'horizontal');
                ylabel('\eta_3', 'FontSize', 11);
                xlim([0, iteration - 1]);
                ylim([0.2 1.2]);
                grid on;
                %
                subplot (4, 1, 4);
                h1 = plot(iteration_arr - 1, x_arr(4, :), 'LineWidth', 2);
                h2 = yline(batch.motor_efficiency(4), '--k', 'Color', 'r', 'LineWidth', 2);
                legend('Estimated', 'True', 'Location', 'southeast', 'Orientation', 'horizontal');
                xlabel('Iterations', 'FontSize', 11);
                ylabel('\eta_4', 'FontSize', 11);
                xlim([0, iteration - 1]);
                ylim([0.2 1.2]);
                grid on;
                exportgraphics(fig, 'efficiency_iteration.png');
                
                fig = figure('Name', 'Trajectory residual norm per Interior-Point Round');
                hold on;
                %
                num_rounds = numel(residual_per_k);
                max_steps  = 0;
                for kk = 1:num_rounds
                    max_steps = max(max_steps, numel(residual_per_k{kk}));
                end
                %
                legends = strings(1, num_rounds);
                for kk = 1:num_rounds
                    rvec = residual_per_k{kk};
                    plot(0:numel(rvec)-1, rvec, 'LineWidth', 1.7);
                    legends(kk) = "k=" + kk;
                end
                %
                xlabel('Interior-point steps', 'FontSize', 11);
                ylabel('$\|\textbf{r(s)}\|$', 'Interpreter', 'latex', 'FontSize', 11);
                xlim([0, max_steps-1]);
                title('Trajectory residual norm per interior-point round', 'FontSize', 11);
                grid on;
                legend(legends, 'Location', 'northeast');
                hold off;
                exportgraphics(fig, 'residual_iteration_multi_rounds.png');
                
                fig = figure('Name', 'KKT dual residual norm vs Iteration');
                % (1) Primal residual
                subplot(3, 1, 1);
                semilogy(iteration_arr(1:end) - 1, residual_arr(1:end), 'o-', ...
                    'LineWidth', 1.7, ...
                    'Color', '#1f77b4', ...
                    'MarkerEdgeColor', '#1f77b4', ...
                    'MarkerFaceColor', '#1f77b4');
                ylabel('$\|\mathbf{r(s)}\|$', 'Interpreter', 'latex', 'FontSize', 11);
                xlim([1, iteration - 1]);
                title('Primal residual norm (log scale)', 'FontSize', 11);
                grid on;
                % (2) Dual residual
                subplot(3, 1, 2);
                semilogy(iteration_arr(1:end) - 1, r_dual_now_arr(1:end), 'o-', ...
                    'LineWidth', 1.7, ...
                    'Color', '#1f77b4', ...
                    'MarkerEdgeColor', '#1f77b4', ...
                    'MarkerFaceColor', '#1f77b4');
                ylabel('$\|\mathbf{r}_{\textbf{dual}}\|$', 'Interpreter', 'latex', 'FontSize', 11);
                xlim([1, iteration - 1]);
                title('Dual residual norm (log scale)', 'FontSize', 11);
                grid on;
                % (3) Duality gap
                subplot(3, 1, 3);
                semilogy(iteration_arr(1:end) - 1, gap_arr(1:end), 'o-', ...
                    'LineWidth', 1.7, ...
                    'Color', '#1f77b4', ...
                    'MarkerEdgeColor', '#1f77b4', ...
                    'MarkerFaceColor', '#1f77b4');
                xlabel('Iterations', 'FontSize', 11);
                ylabel('$\hat{\delta}(s_t,\lambda)$', 'Interpreter', 'latex', 'FontSize', 11);
                xlim([1, iteration - 1]);
                title('Surrogate duality gap (log scale)', 'FontSize', 11);
                grid on;
                exportgraphics(fig, 'kkt_convergence.png');
                
                pause;
                close all;
            end
            
            ret_x = x_cur;
        end
        
        function [ret_x, cond_max] = run(obj, batch, x, x_last)
            [ret_x, cond_max] = obj.run_primal_dual(batch, x, x_last);
        end
    end
end