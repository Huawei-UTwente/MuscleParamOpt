%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct collocation format of the trajectory optimization in muscle dynamics
%
% By: Huawei Wang
% Date: August 2, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [df_dx, df_ddx, df_dp] = ...
    directCollocationDyn_diff_MPO1(x, dx, lmt, par, M, P, S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x: the state variables at current frame
% dx: the derivative state variables at current frame
% lmt: experimental muscle tendon length
% par: parameters of the muscle model (to be optimized)
% M: number of muscles
% P: number of optimizing muscle parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % muscle model parameters

    lce_opt = par(1:M);
    lt_slack = par(M+1:2*M);
    theta0 = par(2*M+1:3*M);

    a = x(1 : M);  % muscle activation at current frame
%     da = x(M+1 : 2*M);  % muscle activation differentiation at current frame
    lce = x(2*M+1 : 3*M);  % length of contract element at current frame
    dlce = x(3*M + 1 : 4*M);  % length of contract element differentiation at current frame

%     d_a = dx(1 : M);  % muscle activation differentiation at current frame
%     d_lce = dx(2*M + 1 : 3*M);  % length of contract element differentiation at current frame
            
    % constraints of finite differetiations
    % fd = [da - d_a, dlce - d_lce];
    
    dfd_dx = zeros(2*M, S*M);
    dfd_ddx = zeros(2*M, S*M);
    
    dfd_dx(1 : M, M+1 : 2*M) = eye(M);
    dfd_ddx(1 : M, 1 : M) = -eye(M);
    
    dfd_dx(M+1 : 2*M, 3*M+1 : 4*M) = eye(M);
    dfd_ddx(M+1 : 2*M, 2*M+1 : 3*M) = -eye(M);

%     % muscle activation dynamics
%     [fa, dfa_da, dfa_dda, dfa_du] ...
%         = activationDyn_Groote_diff(a, da, s_em, M);

    % activiation nonlinearity
    % a_non = activationNonlinearty(a);
    % da_non_da = activationNonlinearty_diff(a);
    a_non = a;
    da_non_da = 1;


    % muscle force dynamics
    [~, dfm_da_non, dfm_dlce, dfm_ddlce, dfm_dlce_opt, dfm_dlt_slack, dfm_dtheta0]...
        = contractionDyn_Groote_diff(lmt, a_non, lce, dlce, lce_opt, lt_slack, theta0);

    % generate the combined constraints
    % f = [fd, fm];
    
    dfm_dx = zeros(M, S*M);
    dfm_dx(:, 1:M) = diag(dfm_da_non.*da_non_da);
    dfm_dx(:, 2*M + 1:3*M) = diag(dfm_dlce);
    dfm_dx(:, 3*M + 1:4*M) = diag(dfm_ddlce);
    
    dfm_dp = zeros(M, P);
    
    dfm_dp(:, 1:M) = diag(dfm_dlce_opt);
    dfm_dp(:, M + 1:2*M) = diag(dfm_dlt_slack);
    dfm_dp(:, 2*M + 1:3*M) = diag(dfm_dtheta0);

    df_dx = [dfd_dx; dfm_dx];
    df_ddx = [dfd_ddx; zeros(M, S*M)];
    df_dp = [zeros(2*M, P); dfm_dp];
            
end