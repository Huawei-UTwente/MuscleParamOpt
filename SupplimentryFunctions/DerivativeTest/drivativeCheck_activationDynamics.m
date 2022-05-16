%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative test of the activation dynamics
%
% By: Huawei Wang
% Date: 12/06/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

% initlize parameters
M = 2;  % number of muscles
a = [0.3, 0.7];  % muscle activation
da = [1.5, 1.3];  % drivative of muscle activation
u = [0.4, 0.6];  % muscle excitation
Ta = [0.015, 0.02];  % activation time constant
Td = [0.060, 0.050];  % deactivation time constant

% calculate activation dynamics derivatives
[~, df_da, df_dda, df_du] =...
    activationDyn_Groote_diff(a, da, u, M);

df_da_equ = diag(df_da);
df_dda_equ = diag(df_dda);
df_du_equ = diag(df_du);
% df_dTa_equ = diag(df_dTa);
% df_dTd_equ = diag(df_dTd);

% finite differentiation
delta = 1e-4;

% differentiation of muscle activation a
for ia = 1:length(a)
    
   delta_a = a(ia)*delta;
   
   % get values with upper change of a
   a(ia) = a(ia) + delta_a;
   f_up = activationDyn_Groote(a, da, u, M);

   % get values with upper change of a
   a(ia) = a(ia) - 2*delta_a;
   f_dw = activationDyn_Groote(a, da, u, M);

    a(ia) = a(ia) + delta_a;

    df_da_fd(:, ia) = (f_up - f_dw)/(2*delta_a);

end

% differentiation of muscle activation drivatives da
for ida = 1:length(da)
    
   delta_da = da(ida)*delta;
   
   % get values with upper change of a
   da(ida) = da(ida)  + delta_da;
   f_up = activationDyn_Groote(a, da, u, M);

   % get values with upper change of a
   da(ida) = da(ida)  - 2*delta_da;
   f_dw = activationDyn_Groote(a, da, u, M);

    da(ida) = da(ida)  + delta_da;

    df_dda_fd(:, ida) = (f_up - f_dw)/(2*delta_da);

end

% differentiation of muscle activation excitation u
for iu = 1:length(u)
    
   delta_u = u(iu)*delta;
    
   % get values with upper change of a
   u(iu) = u(iu) + delta_u;
   f_up = activationDyn_Groote(a, da, u, M);

   % get values with upper change of a
   u(iu) = u(iu) - 2*delta_u;
   f_dw = activationDyn_Groote(a, da, u, M);

    u(iu) = u(iu) + delta_u;

    df_du_fd(:, iu) = (f_up - f_dw)/(2*delta_u);

end

% differentiation of muscle activation time constant Ta
% for iTa = 1:length(Ta)
%     
%    delta_Ta = Ta(iTa)*delta; 
%     
%    % get values with upper change of a
%    Ta(iTa) = Ta(iTa) + delta_Ta;
%    f_up = activationDyn_Groote(a, da, u, M);
% 
%    % get values with upper change of a
%    Ta(iTa) = Ta(iTa) - 2*delta_Ta;
%    f_dw = activationDyn_Groote(a, da, u, M);
% 
%     Ta(iTa) = Ta(iTa) + delta_Ta;
% 
%     df_dTa_fd(:, iTa) = (f_up - f_dw)/(2*delta_Ta);
% 
% end
% 
% % differentiation of muscle deactivation time constant Td
% for iTd = 1:length(Td)
%    
%    delta_Td = Td(iTd)*delta; 
%     
%    % get values with upper change of a
%    Td(iTd) = Td(iTd) + delta_Td;
%    f_up = activationDyn_Groote(a, da, u, M);
% 
%    % get values with upper change of a
%    Td(iTd) = Td(iTd) - 2*delta_Td;
%    f_dw = activationDyn_Groote(a, da, u, M);
% 
%     Td(iTd) = Td(iTd) + delta_Td;
% 
%     df_dTd_fd(:, iTd) = (f_up - f_dw)/(2*delta_Td);
% 
% end

% check the drivative differences
tolerance = 1e-4;

% difference check of df_da
errorid_df_da = diffEvaluate(df_da_equ, df_da_fd, tolerance);
% difference check of df_dda
errorid_df_dda = diffEvaluate(df_dda_equ, df_dda_fd, tolerance);
% difference check of df_du
errorid_df_du = diffEvaluate(df_du_equ, df_du_fd, tolerance);
% difference check of df_dTa
% errorid_df_dTa = diffEvaluate(df_dTa_equ, df_dTa_fd, tolerance);
% % difference check of df_dTd
% errorid_df_dTd = diffEvaluate(df_dTd_equ, df_dTd_fd, tolerance);

if ~isempty(errorid_df_da)
   printf('Differentiations in df_da beyond thresholds\n')
end

if ~isempty(errorid_df_dda)
   printf('Differentiations in df_dda beyond thresholds\n')
end

if ~isempty(errorid_df_du)
   printf('Differentiations in df_du beyond thresholds\n')
end

% if ~isempty(errorid_df_dTa)
%    printf('Differentiations in df_dTa beyond thresholds\n')
% end
% 
% if ~isempty(errorid_df_dTd)
%    printf('Differentiations in df_dTd beyond thresholds\n')
% end




