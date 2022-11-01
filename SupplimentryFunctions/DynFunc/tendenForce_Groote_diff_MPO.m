function [Fse_S, dFse_S_dlce, dFse_S_dlce_opt, dFse_S_dlt_slack, dFse_S_dtheta0]...
            = tendenForce_Groote_diff_MPO(lmt, lce, lce_opt, lt_slack, theta0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This is the tenden force calculation
%
% By: Huawei Wang
% Date: August 1, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % some default parameters, from Groote 2016 paper
    c1 = 0.2;
    c2 = 1;      %0.995;
    c3 = 0.2;    %0.25;
    % kT = 35;   %4 strain
    % kT = 20;   %8.5 strain
%     kT = 25;   %6.5 strain
    kT = 17.5;   %10 strain
    
    % pennation angle calculation and its derivatives
    [cos_theta, dcos_theta_dlce, dcos_theta_dlce_opt, dcos_theta_dtheta0] ...
                = pennationAngSmooth_diff(lce, lce_opt, theta0);
        
    % calcualte tenden lenght
    lt = lmt - lce.*cos_theta;
    
    % calculate derivatives
    dlt_dlce = -cos_theta -lce.*dcos_theta_dlce;
    dlt_dcos_theta = -lce;
    
    % calculate normalized tendon length 
    lt_nor = lt./lt_slack;
    
    dlt_nor_dlt = 1./lt_slack;
    dlt_nor_dlt_slack = -lt./((lt_slack).^2);
    
    % calculate tendon force
    Fse = c1*exp(kT.*(lt_nor - c2)) - c3;
    
    % smooth the Fse & removing negative Fse values
    Fse_S = (sqrt(Fse.^2 + 1e-5) + Fse)/2;
    
    % calculate differentiation
    dFse_dlce = (Fse + c3)*kT.*dlt_nor_dlt.*dlt_dlce;
    dFse_dlce_opt = (Fse + c3)*kT.*dlt_nor_dlt.*dlt_dcos_theta.*dcos_theta_dlce_opt;
    dFse_dlt_slack = (Fse + c3)*kT.*dlt_nor_dlt_slack;
    dFse_dtheta0 = (Fse + c3)*kT.*dlt_nor_dlt.*dlt_dcos_theta.*dcos_theta_dtheta0;
    
    dFse_S_dFse = 0.5*Fse./sqrt(Fse.^2 + 1e-5) + 0.5;
    
    dFse_S_dlce = dFse_S_dFse.*dFse_dlce;
    dFse_S_dlce_opt = dFse_S_dFse.*dFse_dlce_opt;
    dFse_S_dlt_slack = dFse_S_dFse.*dFse_dlt_slack;
    dFse_S_dtheta0 = dFse_S_dFse.*dFse_dtheta0;
end