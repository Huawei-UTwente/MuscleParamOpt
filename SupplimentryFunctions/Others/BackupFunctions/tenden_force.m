function [Fse, dFse_dKs, dFse_dlce, dFse_dlce_opt, dFse_dtheta0,...
   dFse_dlt_slack, dFse_dFmax] = tenden_force(lmt, lce, lce_opt, theta0, lt_slack, Ks, Fmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This is the tenden force calculation
%
% By: Huawei Wang
% Date: August 1, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % angle of the contraction element
    h_ce = lce_opt.*sin(theta0);
    
    % lce^2 - h^2 will always larger than 0 by bounding the lce in
    % optimization. Therefore smoothing is not needed. 
    
    smooth_delta = 1e-8;
    
    y = lce.^2 - h_ce.^2;

    dy_dlce = 2*lce;
    dy_dh = -2*h_ce;

    dy = (sqrt(y.^2 + smooth_delta) + y)/2;
    
    ddy_dy = y./sqrt(y.^2 + smooth_delta)/2 + 1/2; 
    
    ddy_dlce = ddy_dy.*dy_dlce;
    ddy_dh = ddy_dy.*dy_dh;
    
    cos_theta = sqrt(dy)./lce;  % calcualte cos(theta) at current lce length
    
    % calculate its derivatives    
    dcos_theta_ddy = 1./(2*sqrt(dy).*lce);
    
    dcos_theta_dlce = dcos_theta_ddy.*ddy_dlce - sqrt(dy)./lce.^2;
    dcos_theta_dh = dcos_theta_ddy.*ddy_dh;
    
    % calculate derivatives
    dh_dlce_opt = sin(theta0);
    dh_dtheta0 = lce_opt.*cos(theta0);
    
    dcos_theta_dlce_opt = dcos_theta_dh.*dh_dlce_opt;
    dcos_theta_dtheta0 = dcos_theta_dh.*dh_dtheta0;
        
    % calcualte tenden force
    lt = lmt - lce.*cos_theta;
    
    dlt_dlce = -cos_theta -lce.*dcos_theta_dlce;
    dlt_dcos_theta = -lce;
    
    d_lt = lt - lt_slack;
    
    dd_lt = (sqrt(d_lt.^2 + smooth_delta) + d_lt)/2;
    
    ddd_lt_dlt = 1./(2*sqrt(d_lt.^2 + smooth_delta)).*d_lt + 1/2;
    ddd_lt_dlt_slack = -1./(2*sqrt(d_lt.^2 + smooth_delta)).*d_lt - 1/2;
    
    ddd_lt_dlce = ddd_lt_dlt.*dlt_dlce;
    ddd_lt_dcos_theta = ddd_lt_dlt.*dlt_dcos_theta;
    
    
    Fse = 0.1*Fmax.*(lt - lt_slack) + (1/(0.05^2))*Ks.*Fmax.*dd_lt.^2;
    
    % calcualte derivatives
    dFse_dKs = (1/(0.05^2))*Fmax.*dd_lt.^2;
    dFse_dlce = 0.1*Fmax.*dlt_dlce + 2*(1/(0.05^2))*Ks.*Fmax.*dd_lt.*ddd_lt_dlce;
    dFse_dlce_opt = 0.1*Fmax.*dlt_dcos_theta.*dcos_theta_dlce_opt...
        + 2*(1/(0.05^2))*Ks.*Fmax.*dd_lt.*ddd_lt_dcos_theta.*dcos_theta_dlce_opt;
    dFse_dtheta0 = 0.1*Fmax.*dlt_dcos_theta.*dcos_theta_dtheta0...
        + 2*(1/(0.05^2))*Ks.*Fmax.*dd_lt.*ddd_lt_dcos_theta.*dcos_theta_dtheta0;
    dFse_dlt_slack = -0.1*Fmax + 2*(1/(0.05^2))*Ks.*Fmax.*dd_lt.*ddd_lt_dlt_slack;
    dFse_dFmax = 0.1*(lt - lt_slack) + (1/(0.05^2))*Ks.*dd_lt.^2;
   
end