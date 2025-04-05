classdef motor_estimator
    properties
        math;
        
        n = 10; % Number of trajectory points (-1 for # of segments)
        I = eye(3);
        g = 9.8;
        e3 = [0; 0; 1];
        dt;
        dt2; % dt squared
        
        Sigma_i = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0;  % sigma_vx
            0, 1, 0, 0, 0, 0, 0, 0, 0, 0;  % sigma_vy
            0, 0, 1, 0, 0, 0, 0, 0, 0, 0;  % sigma_vz
            0, 0, 0, 1, 0, 0, 0, 0, 0, 0;  % sigma_px
            0, 0, 0, 0, 1, 0, 0, 0, 0, 0;  % sigma_py
            0, 0, 0, 0, 0, 1, 0, 0, 0, 0;  % sigma_pz
            0, 0, 0, 0, 0, 0, 1, 0, 0, 0;  % sigma_Wx
            0, 0, 0, 0, 0, 0, 0, 1, 0, 0;  % sigma_Wy
            0, 0, 0, 0, 0, 0, 0, 0, 1, 0;  % sigma_Wz
            0, 0, 0, 0, 0, 0, 0, 0, 0, 1]; % sigma_R
        Sigma = [];
        Sigma_inv = [];
        
        lower_bound = 0;
        higher_bound = 1.1;
    end
    
    methods
        function ret = init(obj, math, dt)
            obj.math = math;
            obj.dt = dt;
            obj.dt2 = dt^2;
            
            for i = 1:obj.n-1
                obj.Sigma = blkdiag(obj.Sigma, obj.Sigma_i);
            end
            obj.Sigma_inv = inv(obj.Sigma);
            
            ret = obj;
        end
        
        function batch = get_new_batch_random(obj, data)
            idx = randi([1, size(data.time_arr, 2) - obj.n]);
            batch.p = data.pos_arr(:, idx:idx+obj.n);
            batch.v = data.vel_arr(:, idx:idx+obj.n);
            batch.W = data.W_arr(:, idx:idx+obj.n);
            batch.R = data.R_arr(:, :, idx:idx+obj.n);
            batch.f = data.ideal_f_arr(:, idx:idx+obj.n);
            batch.f_motors = data.f_motors_arr(:, idx:idx+obj.n);
            batch.M = data.ideal_M_arr(:, idx:idx+obj.n);
            
            batch.m = data.m;
            batch.J = data.J;
        end
        
        function ret = calc_log_barrier_gradient(obj, x)
            ret = ...
                [1/(obj.higher_bound - x(1)) - 1/(x(1)-obj.lower_bound);
                1/(obj.higher_bound - x(2)) - 1/(x(2)-obj.lower_bound);
                1/(obj.higher_bound - x(3)) - 1/(x(3)-obj.lower_bound);
                1/(obj.higher_bound - x(4)) - 1/(x(4)-obj.lower_bound)];
        end
        
        function ret = calc_log_barrier_hessian(obj, x)
            H11 = 1/((obj.higher_bound + x(1))^2) + 1/((x(1)-obj.lower_bound)^2);
            H22 = 1/((obj.higher_bound + x(2))^2) + 1/((x(2)-obj.lower_bound)^2);
            H33 = 1/((obj.higher_bound + x(3))^2) + 1/((x(3)-obj.lower_bound)^2);
            H44 = 1/((obj.higher_bound + x(4))^2) + 1/((x(4)-obj.lower_bound)^2);
            ret = diag([H11, H22, H33, H44]);
        end
        
        function ret = calc_residual_vector(obj, m, J, x, f_motors, v, p, W, R)
            J_inv = inv(J);
            residual = zeros(10 * (obj.n - 1), 1);
            
            % FIXME
            d = 0.225;
            c = 0.009012;
            
            f = x(1)*f_motors(1) + ...
                x(2)*f_motors(2) + ...
                x(3)*f_motors(3) + ...
                x(4)*f_motors(4);
            M = [d*(-x(1)*f_motors(1) + x(2)*f_motors(2) + x(3)*f_motors(3) - x(4)*f_motors(4)); ...
                d*(+x(1)*f_motors(1) + x(2)*f_motors(2) - x(3)*f_motors(3) - x(4)*f_motors(4)); ...
                c*(-x(1)*f_motors(1) + x(2)*f_motors(2) - x(3)*f_motors(3) + x(4)*f_motors(4))];
            
            for i=1:obj.n-1
                ge3_fRe3m = obj.g*obj.e3 - (f*R(:, :, i)*obj.e3)/m;
                pred_v = v(:, i) + ge3_fRe3m*obj.dt;
                pred_p = p(:, i) + v(:, i)*obj.dt + ge3_fRe3m*obj.dt2/2;
                
                Jinv_M_WJW = J_inv * (M - cross(W(:, i), J*W(:, i)));
                pred_W = W(:, i) + Jinv_M_WJW*obj.dt;
                pred_R = R(:, :, i) * (obj.I + obj.math.hat_map_3x3(W(:, i))*obj.dt);
                
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
        
        function ret = calc_f_jacobian(obj, m, J, f_motors, v, p, W, R)
            X_DIM = 4; % [eta1; eta2; eta3; et4]
            RES_DIM = 10; % [res_v, res_p, res_W, res_R]
            Jf = zeros(RES_DIM*(obj.n - 1), X_DIM);
            
            Jx = J(1,1);
            Jy = J(2,2);
            Jz = J(3,3);
            
            % FIXME
            d = 0.225;
            c = 0.009012;
            
            for i=1:obj.n-1
                delta_R = R(:,:,i).'*R(:,:,i+1);
                r31 = R(3,3,i);
                r32 = R(3,2,i);
                r33 = R(3,3,i);
                
                % d(res_v)/dx
                r_start = 10*(i-1)+1;
                r_end =   10*(i-1)+3;
                Jf(r_start:r_end, 1:4) = obj.dt/m * ...
                    [f_motors(1)*r31, f_motors(2)*r31, f_motors(3)*r31, f_motors(4)*r31; ...
                    f_motors(1)*r32, f_motors(2)*r32, f_motors(3)*r32, f_motors(4)*r32; ...
                    f_motors(1)*r33, f_motors(2)*r33, f_motors(3)*r33, f_motors(4)*r33];
                
                % d(res_p)/dx
                r_start = 10*(i-1)+4;
                r_end =   10*(i-1)+6;
                Jf(r_start:r_end, 1:4) = obj.dt2/(2*m) * ...
                    [f_motors(1)*r31, f_motors(2)*r31, f_motors(3)*r31, f_motors(4)*r31; ...
                    f_motors(1)*r32, f_motors(2)*r32, f_motors(3)*r32, f_motors(4)*r32; ...
                    f_motors(1)*r33, f_motors(2)*r33, f_motors(3)*r33, f_motors(4)*r33];
                
                % d(res_W)/dx
                r_start = 10*(i-1)+7;
                r_end =   10*(i-1)+9;
                Jf(r_start:r_end,1:4) = obj.dt * ...
                    [+f_motors(1)*d/Jx, -f_motors(2)*d/Jx, -f_motors(3)*d/Jx, +f_motors(4)*d/Jx; ...
                    -f_motors(1)*d/Jy, -f_motors(2)*d/Jy, +f_motors(3)*d/Jy, +f_motors(4)*d/Jy; ...
                    +f_motors(1)*c/Jz, -f_motors(2)*c/Jz, +f_motors(3)*c/Jz, -f_motors(4)*c/Jz];
                
                % d(res_R)/dx
                r_start = 10*(i-1)+10;
                r_end =   10*(i-1)+10;
                r1 = -((delta_R(2,3)-delta_R(3,2))*f_motors(1)*d)/Jx ...
                    + ((delta_R(3,1)-delta_R(1,3))*f_motors(1)*d)/Jy ...
                    - ((delta_R(1,2)-delta_R(2,1))*f_motors(1)*c)/Jz;
                r2 = +((delta_R(2,3)-delta_R(3,2))*f_motors(2)*d)/Jx ...
                    + ((delta_R(3,1)-delta_R(1,3))*f_motors(2)*d)/Jy ...
                    + ((delta_R(1,2)-delta_R(2,1))*f_motors(2)*c)/Jz;
                r3 = +((delta_R(2,3)-delta_R(3,2))*f_motors(3)*d)/Jx ...
                    - ((delta_R(3,1)-delta_R(1,3))*f_motors(3)*d)/Jy ...
                    - ((delta_R(1,2)-delta_R(2,1))*f_motors(3)*c)/Jz;
                r4 = -((delta_R(2,3)-delta_R(3,2))*f_motors(4)*d)/Jx ...
                    - ((delta_R(3,1)-delta_R(1,3))*f_motors(4)*d)/Jy ...
                    + ((delta_R(1,2)-delta_R(2,1))*f_motors(4)*c)/Jz;
                Jf(r_start:r_end, 1:4) = obj.dt2 * [r1, r2, r3, r4] / 2;
            end
            
            ret = Jf;
        end
        
        function [ret_x, skip] = gauss_newton_x(obj, iteration, batch, x)
            x0 = x;
            x_last = x;
            skip = 0;
            
            lambda = 0.1;
            mu = 10;
            m = 8; % Constraints numbers  (i.e., lower_bound < x < higher_bound)
            while (m / lambda) > 1e-6 % Outer loop for log barrier control
                while 1 % Inner loop for minimization
                    %fprintf("iteration: %d\n", iteration)
                    
                    % Linearization
                    Jf = lambda * calc_f_jacobian(obj, batch.m, batch.J, batch.f_motors, batch.v, batch.p, batch.W, batch.R);
                    Jf_t = Jf.';
                    
                    % Skip Degenerate Jacobian due to insufficient excitement of the trajectory
                    rcond_JtJ = rcond(Jf.'*obj.Sigma_inv*Jf);
                    if rcond_JtJ < 1e-4
                        fprintf('x: Degenerate Jacobian detected – skip this time step');
                        x = x0;
                        ret_x = x;
                        skip = 1;
                        return;
                    end
                    
                    % Calculate log barrier gradient and hessian
                    g = calc_log_barrier_gradient(obj, x);
                    H = calc_log_barrier_hessian(obj, x);
                    
                    % Calculate Gauss-Newton Step
                    residual_f = lambda * calc_residual_vector(obj, batch.m, batch.J, x, batch.f_motors, batch.v, batch.p, batch.W, batch.R);
                    delta_x = -(Jf_t*obj.Sigma_inv*Jf + H) \ (Jf_t*obj.Sigma_inv*residual_f + g);
                    
                    % Update x
                    x_last = x;
                    x = x + delta_x;
                    
                    %disp(norm(x - x_last));
                    if norm(x - x_last) < 1e-6
                        break;
                    end
                end
                lambda = mu * lambda;
            end
            ret_x = x;
        end
        
        function ret = run(obj, data)
            x_avg = [1; 1; 1; 1];
            alpha = 0.05;
            
            for i = 1:20000
                x = [1; 1; 1; 1];
                
                fprintf("iteration: %d\n", i)
                
                batch = get_new_batch_random(obj, data);
                
                % Run optimization for current trajectory
                [x, skip] = gauss_newton_x(obj, i, batch, x);
                if skip == 1
                    continue;
                end
                
                % Low-pass filter for x
                x_avg = alpha*x + (1-alpha)*x_avg;
                x = x_avg;
                
                disp(x_avg);
                %disp(x)
            end
            
            % Return estimated parameters
            ret = x;
        end
    end
end