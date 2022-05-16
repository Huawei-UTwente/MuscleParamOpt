%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative test of the tendon force
%
% By: Huawei Wang
% Date: 12/06/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

% initlize parameters
M = 2;  % number of muscles
lmt = [0.45, 0.35];  % muscle lengths
lce = [0.08, 0.06];  % drivative of muscle activation
lce_opt = [0.06, 0.07];  % muscle excitation
theta0 = [pi/8, pi/6];  % activation time constant
lt_slack = [0.4, 0.3];  % deactivation time constant

% calculate tendon force derivatives
[~, dFse_dlce, dFse_dlce_opt, dFse_dlt_slack, dFse_dtheta0]...
             = tendenForce_Groote_diff(lmt, lce, lce_opt, lt_slack, theta0);

dFse_dlce_equ = diag(dFse_dlce);
dFse_dlce_opt_equ = diag(dFse_dlce_opt);
dFse_dtheta0_equ = diag(dFse_dtheta0);
dFse_dlt_slack_equ = diag(dFse_dlt_slack);

% finite differentiation
delta = 1e-4;

% differentiation of muscle fiber lengths lce
for ia = 1:length(lce)
    
   % get values with upper change
   delta_lce = lce(ia)*delta;
   
   lce(ia) = lce(ia) + delta_lce;
   Fse_up = tendenForce_Groote(lmt, lce, lce_opt, lt_slack, theta0);
   
   % get values with upper change
   lce(ia) = lce(ia) - 2*delta_lce;
   Fse_dw = tendenForce_Groote(lmt, lce, lce_opt, lt_slack, theta0);

   % change back to the original value
   lce(ia) = lce(ia) + delta_lce;

   % calculate the finite differentiation
   dFse_dlce_fd(:, ia) = (Fse_up - Fse_dw)/(2*delta_lce);

end


% differentiation of muscle optimal fiber lengths lce_opt
for ia = 1:length(lce_opt)
    
   % get values with upper change
   delta_lce_opt = lce_opt(ia)*delta;
   
   lce_opt(ia) = lce_opt(ia) + delta_lce_opt;
   Fse_up = tendenForce_Groote(lmt, lce, lce_opt, lt_slack, theta0);
   
   % get values with upper change
   lce_opt(ia) = lce_opt(ia) - 2*delta_lce_opt;
   Fse_dw = tendenForce_Groote(lmt, lce, lce_opt, lt_slack, theta0);

   % change back to the original value
   lce_opt(ia) = lce_opt(ia) + delta_lce_opt;

   % calculate the finite differentiation
   dFse_dlce_opt_fd(:, ia) = (Fse_up - Fse_dw)/(2*delta_lce_opt);

end


% differentiation of muscle pennation angle theta0
for ia = 1:length(theta0)
    
   % get values with upper change
   delta_theta0 = theta0(ia)*delta;
   
   theta0(ia) = theta0(ia) + delta_theta0;
   Fse_up = tendenForce_Groote(lmt, lce, lce_opt, lt_slack, theta0);
   
   % get values with upper change
   theta0(ia) = theta0(ia) - 2*delta_theta0;
   Fse_dw = tendenForce_Groote(lmt, lce, lce_opt, lt_slack, theta0);

   % change back to the original value
   theta0(ia) = theta0(ia) + delta_theta0;

   % calculate the finite differentiation
   dFse_dtheta0_fd(:, ia) = (Fse_up - Fse_dw)/(2*delta_theta0);

end


% differentiation of muscle tendon slack length lt_slack
for ia = 1:length(lt_slack)
    
   % get values with upper change
   delta_lt_slack = lt_slack(ia)*delta;
   
   lt_slack(ia) = lt_slack(ia) + delta_lt_slack;
   Fse_up = tendenForce_Groote(lmt, lce, lce_opt, lt_slack, theta0);
   
   % get values with upper change
   lt_slack(ia) = lt_slack(ia) - 2*delta_lt_slack;
   Fse_dw = tendenForce_Groote(lmt, lce, lce_opt, lt_slack, theta0);

   % change back to the original value
   lt_slack(ia) = lt_slack(ia) + delta_lt_slack;

   % calculate the finite differentiation
   dFse_dlt_slack_fd(:, ia) = (Fse_up - Fse_dw)/(2*delta_lt_slack);

end

% check the drivative differences
tolerance = 1e-4;

% difference check of dFse_dlce
errorid_dFse_dlce = diffEvaluate(dFse_dlce_equ, dFse_dlce_fd, tolerance);
% difference check of dFse_dlce_opt
errorid_dFse_dlce_opt= diffEvaluate(dFse_dlce_opt_equ, dFse_dlce_opt_fd, tolerance);
% difference check of dFse_dtheta0
errorid_dFse_dtheta0 = diffEvaluate(dFse_dtheta0_equ, dFse_dtheta0_fd, tolerance);
% difference check of dFse_dlt_slack
errorid_dFse_dlt_slack = diffEvaluate(dFse_dlt_slack_equ, dFse_dlt_slack_fd, tolerance);

if ~isempty(errorid_dFse_dlce)
   fprintf('Differentiations in dFse_dlce beyond thresholds\n')
end

if ~isempty(errorid_dFse_dlce_opt)
   fprintf('Differentiations in dFse_dlce_opt beyond thresholds\n')
end

if ~isempty(errorid_dFse_dtheta0)
   fprintf('Differentiations in dFse_dtheta0 beyond thresholds\n')
end

if ~isempty(errorid_dFse_dlt_slack)
   fprintf('Differentiations in dFse_dlt_slack beyond thresholds\n')
end
