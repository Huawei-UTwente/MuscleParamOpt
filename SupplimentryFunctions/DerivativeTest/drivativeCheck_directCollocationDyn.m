%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative test of the direct collocation dynamics
%
% By: Huawei Wang
% Date: 12/06/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

% initlize parameters
T = 2;
L = 1;
M = 4;  % number of muscles
S = 4;
C = 4;
P = 4*M; 

%% generate initial guesses
% initialize the optimizing parameters
mus_a = rand(1, T*L*M);
mus_da = 1 - 2*rand(1, T*L*M);

lce = 0.03 + 0.02*rand(1, T*L*M);
dlce = -0.05 + 0.1*rand(1, T*L*M);

x = [mus_a, mus_da, lce, dlce];

mus_da = 1 - 2*rand(1, T*L*M);
mus_dda = 2 - 4*rand(1, T*L*M);

dlce = -0.05 + 0.1*rand(1, T*L*M);
ddlce = -0.1 + 0.2*rand(1, T*L*M);

dx = [mus_da, mus_dda, dlce, ddlce];

u = rand(1, L*M*T);

lce_opt0 = 0.03 + 0.02*rand(1, M);
lt_slack0 = 0.3 + 0.2*rand(1, M);
theta0 = pi/18 + pi/18*rand(1, M);
Fmax0 = 500 + 500*rand(1, M);

p = [lce_opt0, lt_slack0, theta0, Fmax0];

lmt = 0.45 + 0.15*rand(1, T*L*M);

[f, dfdx, dfddx, dfdu, dfdp] = ...
    directCollocationDyn_diff(x, dx, u, p, T, L, M, S, C, lmt);

% finite differentiation
delta = 1e-6;

% differentiation of input parameters
for ia = 1:length(x)
    
   % get values with upper change
   delta_x = x(ia)*delta;
   
   x(ia) = x(ia) + delta_x;
   f_up = directCollocationDyn(x, dx, u, p, T, L, M, S, C, lmt);
   
   % get values with upper change
   x(ia) = x(ia) - 2*delta_x;
   f_dw = directCollocationDyn(x, dx, u, p, T, L, M, S, C, lmt);

   % change back to the original value
   x(ia) = x(ia) + delta_x;

   % calculate the finite differentiation
   df_dx_fd(:, ia) = (f_up - f_dw)/(2*delta_x);

end

% differentiation of input parameters
for ia = 1:length(dx)
    
   % get values with upper change
   delta_x = dx(ia)*delta;
   
   dx(ia) = dx(ia) + delta_x;
   f_up = directCollocationDyn(x, dx, u, p, T, L, M, S, C, lmt);
   
   % get values with upper change
   dx(ia) = dx(ia) - 2*delta_x;
   f_dw = directCollocationDyn(x, dx, u, p, T, L, M, S, C, lmt);

   % change back to the original value
   dx(ia) = dx(ia) + delta_x;

   % calculate the finite differentiation
   df_ddx_fd(:, ia) = (f_up - f_dw)/(2*delta_x);

end

% differentiation of input parameters
for ia = 1:length(u)
    
   % get values with upper change
   delta_x = x(ia)*delta;
   
   u(ia) = u(ia) + delta_x;
   f_up = directCollocationDyn(x, dx, u, p, T, L, M, S, C, lmt);
   
   % get values with upper change
   u(ia) = u(ia) - 2*delta_x;
   f_dw = directCollocationDyn(x, dx, u, p, T, L, M, S, C, lmt);

   % change back to the original value
   u(ia) = u(ia) + delta_x;

   % calculate the finite differentiation
   df_du_fd(:, ia) = (f_up - f_dw)/(2*delta_x);

end

% differentiation of input parameters
for ia = 1:length(p)
    
   % get values with upper change
   delta_x = p(ia)*delta;
   
   p(ia) = p(ia) + delta_x;
   f_up = directCollocationDyn(x, dx, u, p, T, L, M, S, C, lmt);
   
   % get values with upper change
   p(ia) = p(ia) - 2*delta_x;
   f_dw = directCollocationDyn(x, dx, u, p, T, L, M, S, C, lmt);

   % change back to the original value
   p(ia) = p(ia) + delta_x;

   % calculate the finite differentiation
   df_dp_fd(:, ia) = (f_up - f_dw)/(2*delta_x);

end



% check the drivative differences
tolerance = 1e-4;

% difference check of df_dx
errorid_dfdx = diffEvaluate(dfdx, df_dx_fd, tolerance);
% difference check of df_da
errorid_dfddx = diffEvaluate(dfddx, df_ddx_fd, tolerance);
% difference check of df_da
errorid_dfdu = diffEvaluate(dfdu, df_du_fd, tolerance);
% difference check of df_da
errorid_dfdp = diffEvaluate(dfdp, df_dp_fd, tolerance);

if ~isempty(errorid_dfdx)
   fprintf('Differentiations in dfdx beyond thresholds\n')
end
if ~isempty(errorid_dfddx)
   fprintf('Differentiations in dfddx beyond thresholds\n')
end
if ~isempty(errorid_dfdu)
   fprintf('Differentiations in dfdu beyond thresholds\n')
end
if ~isempty(errorid_dfdp)
   fprintf('Differentiations in dfdp beyond thresholds\n')
end