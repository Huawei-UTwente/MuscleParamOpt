%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct collocation format of the trajectory optimization in muscle dynamics
%
% By: Huawei Wang
% Date: August 2, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [f, dfdx, dfddx, dfdu, dfdp] = ...
    direct_collocation_dyn(x, dx, u, p, T, L, M, S, lmt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x: the state variables at the ith frame
% x1: the state variables at the (i+1)th frame
% u: muscle excitation
% p: parameters of the muscle model (to be optimized)
% T: number of data trails
% L: number of legs
% M: number of muscles
% S: number of states
% lmt: experimental muscle tendon length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    l_row = S*L*M*T;  % number of rows in the sub jacobian matrix
    l_col = length(x);  % number of columns in the sub jacobian matrix
    f = zeros(1, l_row);  % number of elements in the sub constraint vector
    dfdx = zeros(l_row, l_col);  % initialize the sub jacobian matrix dfdx
    dfddx = zeros(l_row, l_col);  % initialize the sub jacobian matrix dfddx
    dfdu = zeros(l_row, M);  % initialize the sub jacobian matrix dfdu
    dfdp = zeros(l_row, length(p));  % initialize the sub jacobian matrix dfdp
    
    % muscle model parameters
    lce_opt = p(1:M);
    lt_slack = p(M+1:2*M);
    theta0 = p(2*M+1:3*M);

    for t = 1:T  % go through all the data trails
        for l = 1:L  % go through all the legs
            
            a = x((t-1)*M*L + (l-1)*M + 1:(t-1)*M*L + l*M);  % muscle activation
            da = x(T*M*L + (t-1)*M*L + (l-1)*M + 1:T*M*L + (t-1)*M*L + l*M);  % muscle activation differentiation
            lce = x(2*T*M*L + (t-1)*M*L + (l-1)*M + 1:2*T*M*L + (t-1)*M*L + l*M);  % length of contract element
            dlce = x(3*T*M*L + (t-1)*M*L + (l-1)*M + 1:3*T*M*L + (t-1)*M*L + l*M);  % length of contract element differentiation

            d_a = dx((t-1)*M*L + (l-1)*M + 1:(t-1)*M*L + l*M);  % muscle activation differentiation
            %dd_a = dx(T*M*L + (t-1)*M*L + (l-1)*M + 1:T*M*L + (t-1)*M*L + l*M);
            d_lce = dx(2*T*M*L + (t-1)*M*L + (l-1)*M + 1:2*T*M*L + (t-1)*M*L + l*M);  % length of contract element differentiation
            %dd_lce = dx(3*T*M*L + (t-1)*M*L + (l-1)*M + 1:3*T*M*L + (t-1)*M*L + l*M);
            
            lmt_ntl = lmt((t-1)*M*L + (l-1)*M + 1:(t-1)*M*L + l*M);  % muscle tendon unit length
            
            u_ntl = u((t-1)*M*L + (l-1)*M + 1:(t-1)*M*L + l*M);  % muscle excitation
            
            % constraints of finite differetiations
            fd = [da - d_a, dlce - d_lce];

            dfd_dda = ones(1, M);
            dfd_dd_a = -ones(1, M);
            dfd_ddlce = ones(1, M);
            dfd_dd_lce = -ones(1, M);

            % muscle activation dynamics
            [fa, dfa_da, dfa_dda, dfa_du] = ...
                activation_dyn_Groote(a, da, u_ntl, M);
            
            % muscle force dynamics
            [fm, dfm_da, dfm_dlce, dfm_ddlce, dfm_dlce_opt, dfm_dlt_slack, dfm_dtheta0]...
                = contraction_dyn_Groote(lmt_ntl, a, lce, dlce, lce_opt, lt_slack, theta0);
         
            % generate the combined constraints
            f((t-1)*M*S*L + (l-1)*M*S + 1:(t-1)*M*S*L + l*M*S) = [fd, fa, fm];
            
            % assign the derivatives of muscle states
            % first 2*num_mus constraints
            dfdx((t-1)*M*L*S + (l-1)*M*S + 1:(t-1)*M*L*S + (l-1)*M*S + M,...
                 T*M*L + (t-1)*M*L + (l-1)*M + 1:T*M*L + (t-1)*M*L + l*M) =...
                 diag(dfd_dda);
            
            dfdx((t-1)*M*L*S + (l-1)*M*S + M + 1:(t-1)*M*L*S + (l-1)*M*S + 2*M,...
                 3*T*M*L + (t-1)*M*L + (l-1)*M + 1:3*T*M*L + (t-1)*M*L + l*M) =...
                 diag(dfd_ddlce);

            dfddx((t-1)*M*L*S + (l-1)*M*S + 1:(t-1)*M*L*S + (l-1)*M*S + M,...
                  (t-1)*M*L + (l-1)*M + 1:(t-1)*M*L + l*M) = diag(dfd_dd_a);
              
            dfddx((t-1)*M*L*S + (l-1)*M*S + M + 1:(t-1)*M*L*S + (l-1)*M*S + 2*M,...
                  2*T*M*L + (t-1)*M*L + (l-1)*M + 1:2*T*M*L + (t-1)*M*L + l*M)=...
                  diag(dfd_dd_lce);

            % from the 2*num_mus to the 3*num_mus constraints
            dfdx((t-1)*M*L*S + (l-1)*M*S + 2*M + 1:(t-1)*M*L*S + (l-1)*M*S + 3*M,...
                 (t-1)*M*L + (l-1)*M + 1:(t-1)*M*L + l*M) = diag(dfa_da);
             
            dfdx((t-1)*M*L*S + (l-1)*M*S + 2*M + 1:(t-1)*M*L*S + (l-1)*M*S + 3*M,...
                 T*M*L + (t-1)*M*L + (l-1)*M + 1:T*M*L + (t-1)*M*L + l*M) =...
                 diag(dfa_dda);

            dfdu((t-1)*M*L*S + (l-1)*M*S + 2*M + 1:(t-1)*M*L*S + (l-1)*M*S + 3*M,...
                 (t-1)*M*L + (l-1)*M + 1:(t-1)*M*L + l*M) = diag(dfa_du);

            dfdp((t-1)*M*L*S + (l-1)*M*S + 2*M + 1:(t-1)*M*L*S + (l-1)*M*S + 3*M,...
                 1) = dfa_dC1;
            dfdp((t-1)*M*L*S + (l-1)*M*S + 2*M + 1:(t-1)*M*L*S + (l-1)*M*S + 3*M,...
                 2) = dfa_dC2;

            % from the 3*num_mus to the 4*num_mus constraints
            dfdx((t-1)*M*L*S + (l-1)*M*S + 3*M + 1:(t-1)*M*L*S + (l-1)*M*S + 4*M,...
                 (t-1)*M*L + (l-1)*M + 1:(t-1)*M*L + l*M) = diag(dfm_da);
             
            dfdx((t-1)*M*L*S + (l-1)*M*S + 3*M + 1:(t-1)*M*L*S + (l-1)*M*S + 4*M,...
                 2*T*M*L + (t-1)*M*L + (l-1)*M + 1:2*T*M*L + (t-1)*M*L + l*M)...
                 = diag(dfm_dlce);
            dfdx((t-1)*M*L*S + (l-1)*M*S + 3*M + 1:(t-1)*M*L*S + (l-1)*M*S + 4*M,...
                 3*T*M*L + (t-1)*M*L + (l-1)*M + 1:3*T*M*L + (t-1)*M*L + l*M)...
                 = diag(dfm_ddlce);

            % derivatives of muslce parameters
            dfdp((t-1)*M*L*S + (l-1)*M*S + 3*M +1:(t-1)*M*L*S + (l-1)*M*S + 4*M,...
                1:M) = diag(dfm_dlce_opt);
            dfdp((t-1)*M*L*S + (l-1)*M*S + 3*M +1:(t-1)*M*L*S + (l-1)*M*S + 4*M,...
                M+1:2*M) = diag(dfm_dlt_slack);
            dfdp((t-1)*M*L*S + (l-1)*M*S + 3*M +1:(t-1)*M*L*S + (l-1)*M*S + 4*M,...
                2*M+1:3*M) = diag(dfm_dtheta0);
            
        end
    end
end