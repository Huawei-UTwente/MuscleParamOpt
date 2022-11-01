%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative test of the contaction dynamics
%
% By: Huawei Wang
% Date: 12/06/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

% initlize parameters
M = 2;  % number of muscles
lmt = [0.45, 0.35];  % muscle lengths
a = [0.4, 0.7];      % muscle activation
lce = [0.08, 0.06];  % muscle fiber length
dlce = [-0.4, 0.3];  % derivative of muscle fiber length
lce_opt = [0.06, 0.07];  % muscle optimal fiber length
theta0 = [pi/8, pi/6];  % pennation angle
lt_slack = [0.4, 0.3];  % tendon slack length

% calculate tendon force derivatives
[~, df_da, df_dlce, df_ddlce, df_dlce_opt, df_dlt_slack, df_dtheta0]...
         = contractionDyn_Groote_S_diff(lmt, a, lce, dlce, lce_opt, lt_slack, theta0);

df_da_equ = diag(df_da);
df_dlce_equ = diag(df_dlce);
df_ddlce_equ = diag(df_ddlce);
df_dlce_opt_equ = diag(df_dlce_opt);
df_dlt_slack_equ = diag(df_dlt_slack);
df_dtheta0_equ = diag(df_dtheta0);

% finite differentiation
delta = 1e-4;

% differentiation of muscle fiber lengths lce
for ia = 1:length(a)
    
   % get values with upper change
   delta_a = a(ia)*delta;
   
   a(ia) = a(ia) + delta_a;
   f_up = contractionDyn_Groote_S(lmt, a, lce, dlce, lce_opt, lt_slack, theta0);
   
   % get values with upper change
   a(ia) = a(ia) - 2*delta_a;
   f_dw = contractionDyn_Groote_S(lmt, a, lce, dlce, lce_opt, lt_slack, theta0);

   % change back to the original value
   a(ia) = a(ia) + delta_a;

   % calculate the finite differentiation
   df_da_fd(:, ia) = (f_up - f_dw)/(2*delta_a);

end


% differentiation of muscle fiber lengths lce
for ia = 1:length(a)
    
   % get values with upper change
   delta_a = a(ia)*delta;
   
   a(ia) = a(ia) + delta_a;
   f_up = contractionDyn_Groote_S(lmt, a, lce, dlce, lce_opt, lt_slack, theta0);
   
   % get values with upper change
   a(ia) = a(ia) - 2*delta_a;
   f_dw = contractionDyn_Groote_S(lmt, a, lce, dlce, lce_opt, lt_slack, theta0);

   % change back to the original value
   a(ia) = a(ia) + delta_a;

   % calculate the finite differentiation
   df_da_fd(:, ia) = (f_up - f_dw)/(2*delta_a);

end


% differentiation of muscle fiber lengths lce
for ia = 1:length(lce)
    
   % get values with upper change
   delta_a = lce(ia)*delta;
   
   lce(ia) = lce(ia) + delta_a;
   f_up = contractionDyn_Groote_S(lmt, a, lce, dlce, lce_opt, lt_slack, theta0);
   
   % get values with upper change
   lce(ia) = lce(ia) - 2*delta_a;
   f_dw = contractionDyn_Groote_S(lmt, a, lce, dlce, lce_opt, lt_slack, theta0);

   % change back to the original value
   lce(ia) = lce(ia) + delta_a;

   % calculate the finite differentiation
   df_dlce_fd(:, ia) = (f_up - f_dw)/(2*delta_a);

end


% differentiation of muscle fiber lengths lce
for ia = 1:length(dlce)
    
   % get values with upper change
   delta_a = dlce(ia)*delta;
   
   dlce(ia) = dlce(ia) + delta_a;
   f_up = contractionDyn_Groote_S(lmt, a, lce, dlce, lce_opt, lt_slack, theta0);
   
   % get values with upper change
   dlce(ia) = dlce(ia) - 2*delta_a;
   f_dw = contractionDyn_Groote_S(lmt, a, lce, dlce, lce_opt, lt_slack, theta0);

   % change back to the original value
   dlce(ia) = dlce(ia) + delta_a;

   % calculate the finite differentiation
   df_ddlce_fd(:, ia) = (f_up - f_dw)/(2*delta_a);

end


% differentiation of muscle fiber lengths lce
for ia = 1:length(lce_opt)
    
   % get values with upper change
   delta_a = lce_opt(ia)*delta;
   
   lce_opt(ia) = lce_opt(ia) + delta_a;
   f_up = contractionDyn_Groote_S(lmt, a, lce, dlce, lce_opt, lt_slack, theta0);
   
   % get values with upper change
   lce_opt(ia) = lce_opt(ia) - 2*delta_a;
   f_dw = contractionDyn_Groote_S(lmt, a, lce, dlce, lce_opt, lt_slack, theta0);

   % change back to the original value
   lce_opt(ia) = lce_opt(ia) + delta_a;

   % calculate the finite differentiation
   df_dlce_opt_fd(:, ia) = (f_up - f_dw)/(2*delta_a);

end


% differentiation of muscle fiber lengths lce
for ia = 1:length(lt_slack)
    
   % get values with upper change
   delta_a = lt_slack(ia)*delta;
   
   lt_slack(ia) = lt_slack(ia) + delta_a;
   f_up = contractionDyn_Groote_S(lmt, a, lce, dlce, lce_opt, lt_slack, theta0);
   
   % get values with upper change
   lt_slack(ia) = lt_slack(ia) - 2*delta_a;
   f_dw = contractionDyn_Groote_S(lmt, a, lce, dlce, lce_opt, lt_slack, theta0);

   % change back to the original value
   lt_slack(ia) = lt_slack(ia) + delta_a;

   % calculate the finite differentiation
   df_dlt_slack_fd(:, ia) = (f_up - f_dw)/(2*delta_a);

end


% differentiation of muscle fiber lengths lce
for ia = 1:length(theta0)
    
   % get values with upper change
   delta_a = theta0(ia)*delta;
   
   theta0(ia) = theta0(ia) + delta_a;
   f_up = contractionDyn_Groote_S(lmt, a, lce, dlce, lce_opt, lt_slack, theta0);
   
   % get values with upper change
   theta0(ia) = theta0(ia) - 2*delta_a;
   f_dw = contractionDyn_Groote_S(lmt, a, lce, dlce, lce_opt, lt_slack, theta0);

   % change back to the original value
   theta0(ia) = theta0(ia) + delta_a;

   % calculate the finite differentiation
   df_dtheta0_fd(:, ia) = (f_up - f_dw)/(2*delta_a);

end

% check the drivative differences
tolerance = 1e-4;

% difference check of df_da
errorid_df_da = diffEvaluate(df_da_equ, df_da_fd, tolerance);
% difference check of dFse_dlce
errorid_df_dlce = diffEvaluate(df_dlce_equ, df_dlce_fd, tolerance);
% difference check of dFse_dlce
errorid_df_ddlce = diffEvaluate(df_ddlce_equ, df_ddlce_fd, tolerance);
% difference check of dFse_dlce_opt
errorid_df_dlce_opt = diffEvaluate(df_dlce_opt_equ, df_dlce_opt_fd, tolerance);
% difference check of dFse_dlt_slack
errorid_df_dlt_slack = diffEvaluate(df_dlt_slack_equ, df_dlt_slack_fd, tolerance);
% difference check of dFse_dtheta0
errorid_df_dtheta0 = diffEvaluate(df_dtheta0_equ, df_dtheta0_fd, tolerance);

if ~isempty(errorid_df_da)
   fprintf('Differentiations in df_da beyond thresholds\n')
end

if ~isempty(errorid_df_dlce)
   fprintf('Differentiations in df_dlce beyond thresholds\n')
end

if ~isempty(errorid_df_ddlce)
   fprintf('Differentiations in df_ddlce beyond thresholds\n')
end

if ~isempty(errorid_df_dlce_opt)
   fprintf('Differentiations in df_dlce_opt beyond thresholds\n')
end

if ~isempty(errorid_df_dlt_slack)
   fprintf('Differentiations in df_dlt_slack beyond thresholds\n')
end

if ~isempty(errorid_df_dtheta0)
   fprintf('Differentiations in dFse_df_dtheta0 beyond thresholds\n')
end
