classdef motor_estimator_ekf < handle
    properties
        math
        
        n = 2
        dt;
        dt2;
        
        % constants / params
        g  = 9.8;
        e3 = [0;0;1];
        I3 = eye(3);
        
        % vehicle params
        m;
        J;
        J_inv;
        d;
        c;
        
        % Noises
        Q;          % 22x22 process noise
        R_meas;     % 10x10 measurement noise: [p(3) v(3) W(3) phi(1)]
        
        % State (x = [p(3); v(3); W(3); vec(R)(9); eta(4)])
        x;          % 22x1
        P;          % 22x22
        
        inited = false;
        cond_threshold = 1e4;
    end
    
    methods
        function obj = init(obj, math, dt)
            obj.math = math;
            
            obj.dt  = dt;
            obj.dt2 = dt^2;
            
            % Process noise
            Qp   = 1e-6 * eye(3);
            Qv   = 1e-3 * eye(3);
            Qw   = 1e-4 * eye(3);
            QR   = 1e-5 * eye(9);
            Qeta = 5e-1 * eye(4);
            obj.Q = blkdiag(Qp,Qv,Qw,QR,Qeta); % 22x22
            
            % Measurement noise
            Rp   = 5e-2 * eye(3);
            Rv   = 5e-2 * eye(3);
            Rw   = 5e-2 * eye(3);
            Rphi = 5e-2;
            obj.R_meas = blkdiag(Rp, Rv, Rw, Rphi); % 10x10
            
            % Process covariance matrix
            obj.P = 1e-2 * eye(22);
        end
        
        function [ret_eta, cond] = run(obj, batch, x_init, x_last)
            if ~obj.inited
                obj.m = batch.m;
                obj.J = batch.J;
                obj.J_inv = diag([1/batch.J(1,1), 1/batch.J(2,2), 1/batch.J(3,3)]);
                obj.d = batch.d;
                obj.c = batch.c;
                
                p0 = batch.p(:,1);
                v0 = batch.v(:,1);
                W0 = batch.W(:,1);
                R0 = batch.R(:,:,1);
                eta0 = x_last(:);
                obj.x = obj.pack_x(p0, v0, W0, R0, eta0);
                obj.inited = true;
            end
            
            % Prediction
            fmot_k = batch.f_motors(:,1);
            obj.predict(fmot_k);
            
            % Update
            meas.p = batch.p(:,2);
            meas.v = batch.v(:,2);
            meas.W = batch.W(:,2);
            meas.R = batch.R(:,:,2);
            [cond, ~] = obj.update(meas);
            
            ret_eta = obj.x(19:22);
        end
        
        function predict(obj, fmot)
            [p,v,W,R,eta] = obj.unpack_x(obj.x);
            dt  = obj.dt; dt2 = obj.dt2;
            
            % Velocity and position prediction
            f = eta.'*fmot;
            a = obj.g*obj.e3 - (f/obj.m)*(R*obj.e3);
            v_pred = v + a*dt;
            p_pred = p + v*dt + 0.5*a*dt2;
            
            % Angular velocity prediction
            M = [ obj.d*(-eta(1)*fmot(1) + eta(2)*fmot(2) + eta(3)*fmot(3) - eta(4)*fmot(4));
                obj.d*( eta(1)*fmot(1) + eta(2)*fmot(2) - eta(3)*fmot(3) - eta(4)*fmot(4));
                obj.c*(-eta(1)*fmot(1) + eta(2)*fmot(2) - eta(3)*fmot(3) + eta(4)*fmot(4)) ];
            Wdot = obj.J_inv * (M - cross(W, obj.J*W));
            W_pred = W + Wdot*dt;
            
            % Rotation matrix prediction with first-order update
            R_pred = R + R*obj.math.hat_map_3x3(W)*dt;
            % Project rotation matrix back to SO(3)
            R_pred = obj.project_SO3(R_pred);
            
            % eta prediction
            eta_pred = eta;
            
            % Assemble full-state
            x_pred = obj.pack_x(p_pred, v_pred, W_pred, R_pred, eta_pred);
            
            % Build state transition Jacobian F
            F = eye(22);
            
            % E: 3rd column of vec(R)
            E = zeros(3,9);
            E(:,7:9) = eye(3);
            
            % dv/dR
            F(4:6, 10:18) = -(f/obj.m) * E * dt;
            % dp/dR
            F(1:3, 10:18) = -(f/obj.m) * E * (dt2/2);
            
            % dv/deta
            F(4:6, 19:22) = -(R*obj.e3/obj.m) * (fmot(:)') * dt;
            
            % dp/deta
            F(1:3, 19:22) = -(R*obj.e3/obj.m) * (fmot(:)') * (dt2/2);
            
            % dW/deta
            dM_deta = [-obj.d*fmot(1), +obj.d*fmot(2), +obj.d*fmot(3), -obj.d*fmot(4);
                +obj.d*fmot(1), +obj.d*fmot(2), -obj.d*fmot(3), -obj.d*fmot(4);
                -obj.c*fmot(1), +obj.c*fmot(2), -obj.c*fmot(3), +obj.c*fmot(4) ];
            F(7:9, 19:22) = obj.J_inv * dM_deta * dt;
            
            % F_RR = I_9 + (hat(W)^T ⊗ I_3) dt
            HWT = obj.math.hat_map_3x3(W).';
            F_RR = eye(9) + kron(HWT, eye(3)) * dt;
            F(10:18, 10:18) = F_RR;
            
            % F_R_Wi = vec(R * E_i * dt)
            Ewx = [0 0 0; 0 0 -1; 0 1 0];
            Ewy = [0 0 1; 0 0 0; -1 0 0];
            Ewz = [0 -1 0; 1 0 0; 0 0 0];
            F(10:18, 7) = reshape(R*Ewx*dt, 9, 1);
            F(10:18, 8) = reshape(R*Ewy*dt, 9, 1);
            F(10:18, 9) = reshape(R*Ewz*dt, 9, 1);
            
            % Process covariance predict
            obj.P = F*obj.P*F.' + obj.Q;
            
            % Write back
            obj.x = x_pred;
        end
        
        function [condS, S] = update(obj, meas)
            [p, v, W, R, ~] = obj.unpack_x(obj.x);
            
            z=[];
            h=[];
            H=[];
            Rblk=[];
            
            % Position
            z = [z; meas.p];
            h = [h; p];
            Hp = zeros(3,22);
            Hp(:,1:3) = eye(3);
            H = [H; Hp];
            Rblk = blkdiag(Rblk, obj.R_meas(1:3, 1:3));
            
            % Velocity
            z = [z; meas.v];
            h = [h; v];
            Hv = zeros(3, 22);
            Hv(:, 4:6) = eye(3);
            H = [H; Hv];
            Rblk = blkdiag(Rblk, obj.R_meas(4:6, 4:6));
            
            % Angular velocity
            z = [z; meas.W];
            h = [h; W];
            Hw = zeros(3,22);
            Hw(:,7:9) = eye(3);
            H = [H; Hw];
            Rblk = blkdiag(Rblk, obj.R_meas(7:9, 7:9));
            
            % Attitude (scalar, trace-based angle)
            Rm = meas.R;
            c  = max(-1, min(1, (trace(R.'*Rm)-1)/2));
            phi = acos(c);
            
            z = [z; phi];
            h = [h; 0];
            
            % d phi / d vec(R) = - vec(Rm)^T / (2*sqrt(1-c^2))
            den = 2*sqrt(max(1e-12, 1 - c^2));
            dphi_dR = -(Rm(:)).' / den; % 1x9
            Hphi = zeros(1, 22);
            Hphi(1,10:18) = dphi_dR;
            H = [H; Hphi];
            
            Rphi = obj.R_meas(10,10);
            Rblk = blkdiag(Rblk, Rphi);
            
            % EKF update
            S = H*obj.P*H.' + Rblk;
            
            % Calculate conditional number
            rcondS = rcond(S);
            condS = 1 / max(rcondS, eps);
            
            K = obj.P*H.' / (S + 1e-12*eye(size(S)));
            y = z - h;
            dx = K*y;
            
            % State add
            obj.x = obj.x + dx;
            
            % Reshape and project R back to SO(3)
            [p2, v2, W2, R2, eta2] = obj.unpack_x(obj.x);
            R2 = reshape(R2, 3, 3);
            R2 = obj.project_SO3(R2);
            
            % clamp eta to [1e-3, 1-1e-3]
            eta2 = min(max(eta2, 1e-3), 1-1e-3);
            obj.x = obj.pack_x(p2, v2, W2, R2, eta2);
            
            % covariance
            I = eye(22);
            obj.P = (I - K*H) * obj.P;
            obj.P = (obj.P + obj.P.')/2; % numerical symmetry fix
        end
        
        function [p, v, W, R, eta] = unpack_x(~, x)
            p = x(1:3);
            v = x(4:6);
            W = x(7:9);
            R = reshape(x(10:18), 3, 3);
            eta = x(19:22);
        end
        
        function x = pack_x(~, p, v, W, R, eta)
            x = [p; v; W; R(:); eta];
        end
        
        function R = project_SO3(~, R)
            [U, ~, V] = svd(R);
            R = U * diag([1,1,det(U*V')]) * V';
        end
    end
end