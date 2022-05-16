%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct collocation format of the trajectory optimization in muscle dynamics
%
% By: Huawei Wang
% Date: August 2, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function f = directCollocationDyn_MPO2(x, dx, x_em, lmt, p, M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x: the state variables at current frame
% dx: the derivative state variables at current frame
% x_em: the state variables at the frame that caused by stimulation delay
% lmt: experimental muscle tendon length
% p: parameters of the muscle model (to be optimized)
% M: number of muscles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % muscle model parameters
    lce_opt = p(1:M);
    lt_slack = p(M+1:2*M);
    theta0 = p(2*M+1:3*M);

    a = x(1 : M);  % muscle activation at current frame
    da = x(M+1 : 2*M);  % muscle activation differentiation at current frame
    lce = x(2*M+1 : 3*M);  % length of contract element at current frame
    dlce = x(3*M + 1 : 4*M);  % length of contract element differentiation at current frame

    d_a = dx(1 : M);  % muscle activation differentiation at current frame
    d_lce = dx(2*M + 1 : 3*M);  % length of contract element differentiation at current frame
    
    s_em = x_em(4*M+1 : 5*M);  % the optimizing free muscle stimulation

    % constraints of finite differetiations
    fd = [da - d_a, dlce - d_lce];

    % muscle activation dynamics
    fa = activationDyn_Groote(a, da, s_em);

    % activiation nonlinearity, [not needed for the leg muscles]
    % a_non = activationNonlinearty(a);
    a_non = a;
    
    % muscle force dynamics
    fm = contractionDyn_Groote(lmt, a_non, lce, dlce, lce_opt, ...
                               lt_slack, theta0);

    % generate the combined constraints
    f = [fd, fa, fm];
end