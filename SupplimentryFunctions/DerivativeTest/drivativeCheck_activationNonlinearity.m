%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative check of the activation nonlinearty


clc
clear
close all

T = 2;
M = 4;

u = rand(1, T*M);

delta = 1e-6;

da_du = activationNonlinearty_diff(u);

for i = 1:length(u)
    
    delta_u = u(i)*delta;
    
    u(i) = u(i) + delta_u;
    f_up = activationNonlinearty(u);
    
    u(i) = u(i) - 2*delta_u;
    f_dw = activationNonlinearty(u);
    
    u(i) = u(i) + delta_u;
    
    da_du_fd(:, i) = (f_up - f_dw)/(2*delta_u);
    
end

% check the drivative differences
tolerance = 1e-4;

% difference check of df_da
errorid_da_du = diffEvaluate(da_du, diag(da_du_fd)', tolerance);

if ~isempty(errorid_da_du)
   fprintf('Differentiations in errorid_da_du beyond thresholds\n')
end



