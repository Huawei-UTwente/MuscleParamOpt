
function [f, df_dlce, df_dlce_opt, df_dtheta0, df_dW, df_ddlce,...
     df_dvmax_nor, df_dgmax, df_dA, df_dlce_slack_nor, df_dKp, df_da,...
     df_dKs, df_dlt_slack] = contraction_dyn1_smooth_1(a, lmt, lce,...
     dlce, lce_opt, vmax_nor, gmax, W, A, Kp, lce_slack_nor, Ks,...
     lt_slack, theta0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the contraction dynamics of Hill's muscle model
%
% By: Huawei Wang
% Date: August 1, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fce = exp(-((lce - lce_opt)./(W.*lce_opt)).^2);  % force-length relationship
    
    % calculate derivatives
    dfce_dlce = -2*fce.*(lce - lce_opt)./((W.*lce_opt).^2);
    dfce_dlce_opt = 2*fce.*(lce - lce_opt).*lce./(W.^2.*lce_opt.^3);
    dfce_dW = 2*fce.*(lce - lce_opt).^2./(W.^3.*lce_opt.^2);
    
    Nm = length(dlce);
    
    gdce = zeros(1, Nm);
    dgdce_dvmax = zeros(1, Nm);
    dgdce_ddlce = zeros(1, Nm);
    dgdce_dA = zeros(1, Nm);
    dgdce_dgmax = zeros(1, Nm);
    
    vmax = vmax_nor.*lce_opt;
    dvmax_dvmax_nor = lce_opt;
    dvmax_dlce_opt = vmax_nor;
    
    for m = 1:Nm
        if dlce(m) <= 0
            %force-velocity relationship (Hill-Katz model)
            gdce(m) = (vmax(m) + dlce(m))/(vmax(m) - dlce(m)/A(m));
            % calculate derivatives
            dgdce_dvmax(m) = -(dlce(m)/A(m) + dlce(m))/(vmax(m) - dlce(m)/A(m))^2;
            dgdce_ddlce(m) = (vmax(m)/A(m) + vmax(m))/(vmax(m) - dlce(m)/A(m))^2;
            dgdce_dA(m) = -dlce(m)/A(m)^2*(vmax(m) + dlce(m))/(vmax(m) - dlce(m)/A(m))^2;
            
        else
            c = vmax(m)*A(m)*(gmax(m) - 1)/(A(m) + 1);
            gdce(m) = (gmax(m)*dlce(m) + c)/(dlce(m) + c); % force-velocity relationship (Hill-Katz model)
            
            % calculate derivatives
            dc_dvmax = A(m)*(gmax(m) - 1)/(A(m) + 1);
            dc_dgmax = vmax(m)*A(m)/(A(m) + 1);
            dc_dA = vmax(m)*(gmax(m) - 1)/(A(m) + 1)^2;
            
            dgdce_dc = (dlce(m) - gmax(m)*dlce(m))/(dlce(m) + c)^2;
            
            dgdce_dgmax(m) = dlce(m)/(dlce(m) + c) + dgdce_dc*dc_dgmax;
            dgdce_ddlce(m) = (gmax(m) - 1)*c/(dlce(m) + c)^2;
            dgdce_dvmax(m) = dgdce_dc*dc_dvmax;
            dgdce_dA(m) = dgdce_dc*dc_dA;
        end
    end
    
    dgdce_dvmax_nor = dgdce_dvmax.*dvmax_dvmax_nor;
    dgdce_dlce_opt = dgdce_dvmax.*dvmax_dlce_opt;
        
    % passive contraction element force, delta_lce is to make the force continues
    
    smooth_delta = 1e-8;
    
    lce_slack = lce_slack_nor.*lce_opt;
    
    dlce_slack_dlce_slack_nor = lce_opt;
    dlce_slack_dlce_opt = lce_slack_nor;
    
    d_lce = lce - lce_slack;
    
    dd_lce = (sqrt(d_lce.^2 + smooth_delta) + d_lce)/2;
    
    ddd_lce_dlce = 1./(2*sqrt(d_lce.^2 + smooth_delta)).*d_lce + 1/2;
    ddd_lce_dlce_slack = -1./(2*sqrt(d_lce.^2 + smooth_delta)).*d_lce - 1/2;
    
    ddd_lce_dlce_slack_nor = ddd_lce_dlce_slack.*dlce_slack_dlce_slack_nor;
    ddd_lce_dlce_opt = ddd_lce_dlce_slack.*dlce_slack_dlce_opt;
        
    Fpee = 0.1*(lce - lce_slack) + Kp.*dd_lce.^2./lce_opt.^2;
    
    dFpee_dlce = 0.1 + 2*Kp.*dd_lce./lce_opt.^2.*ddd_lce_dlce;
    dFpee_dlce_slack_nor = -0.1*dlce_slack_dlce_slack_nor + 2*Kp.*dd_lce./lce_opt.^2.*ddd_lce_dlce_slack_nor;
    dFpee_dlce_opt = -0.1*dlce_slack_dlce_opt + 2*Kp.*dd_lce./lce_opt.^2.*ddd_lce_dlce_opt ...
                    - 2*Kp.*dd_lce.^2./lce_opt.^3;
    dFpee_dKp = dd_lce.^2./lce_opt.^2;
    
    % angle of the contraction element
    h = lce_opt.*sin(theta0);
    
    % calcualte the ce length in the parallel direction
    y = lce.^2 - h.^2; 
    
    % calcualte it's derivatives
    dy_dlce = 2*lce;
    dy_dh = -2*h;
    
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
    
    % calcualte the force of the contraction element and PEE together
    Fce = (a.*fce.*gdce + dlce./(1000*lce_opt) + Fpee).*cos_theta;
    
    % calculate overall derivatives
    
    dFce_dlce = (a.*gdce.*dfce_dlce + dFpee_dlce).*cos_theta... 
                 + (a.*fce.*gdce + dlce./(1000*lce_opt) + Fpee).*dcos_theta_dlce;
                 
    dFce_dlce_opt = (a.*dfce_dlce_opt.*gdce + a.*fce.*dgdce_dlce_opt ...
                   - dlce./(1000*lce_opt.^2) + dFpee_dlce_opt).*cos_theta... 
                   + (a.*fce.*gdce + dlce./(1000*lce_opt) + Fpee).*dcos_theta_dlce_opt;
    
    dFce_dW = a.*dfce_dW.*gdce.*cos_theta;
    
    dFce_ddlce = (a.*fce.*dgdce_ddlce + 1./(1000*lce_opt)).*cos_theta;
    
    dFce_dvmax_nor = (a.*fce.*dgdce_dvmax_nor).*cos_theta;
    
    dFce_dgmax = a.*fce.*dgdce_dgmax.*cos_theta;
    
    dFce_dA = a.*fce.*dgdce_dA.*cos_theta;
    
    dFce_dlce_slack_nor = dFpee_dlce_slack_nor.*cos_theta;
    
    dFce_dKp = dFpee_dKp.*cos_theta;
    
    dFce_dtheta0 = (a.*fce.*gdce + dlce./(1000*lce_opt) + Fpee).*dcos_theta_dtheta0;
    
    dFce_da = fce.*gdce.*cos_theta;
        
    % calcualte tenden force
    lt = lmt - lce.*cos_theta;
    
    dlt_dlce = -cos_theta - lce.*dcos_theta_dlce;
    dlt_dcos_theta = -lce;
    
    % smooth the deta lt respect to lt_slack
    d_lt = lt - lt_slack;
    
    dd_lt = (sqrt(d_lt.^2 + smooth_delta) + d_lt)/2;
    
    ddd_lt_dlt = 1./(2*sqrt(d_lt.^2 + smooth_delta)).*d_lt + 1/2;
    ddd_lt_dlt_slack = -1./(2*sqrt(d_lt.^2 + smooth_delta)).*d_lt - 1/2;
    
    ddd_lt_dlce = ddd_lt_dlt.*dlt_dlce;
    ddd_lt_dcos_theta = ddd_lt_dlt.*dlt_dcos_theta;
    
    Fse = 0.1*(lt - lt_slack) + (1/(0.05^2))*Ks.*dd_lt.^2;
    
    % calcualte derivatives
    dFse_dKs = (1/(0.05^2))*dd_lt.^2;
    dFse_dlce = 0.1*dlt_dlce + 2*(1/(0.05^2))*Ks.*dd_lt.*ddd_lt_dlce;
    dFse_dlce_opt = 0.1*dlt_dcos_theta.*dcos_theta_dlce_opt...
        + 2*(1/(0.05^2))*Ks.*dd_lt.*ddd_lt_dcos_theta.*dcos_theta_dlce_opt;
    dFse_dtheta0 = 0.1*dlt_dcos_theta.*dcos_theta_dtheta0...
        + 2*(1/(0.05^2))*Ks.*dd_lt.*ddd_lt_dcos_theta.*dcos_theta_dtheta0;
    
    % dFse_dlmt = Ks*2.*delta_lt.*ddelta_lt_dlt;   % does not need now
    dFse_dlt_slack = -0.1 + 2*(1/(0.05^2))*Ks.*dd_lt.*ddd_lt_dlt_slack;
    
    % force of the contraction element should equal to the force of the tense unit
    f = (Fce - Fse);
    
    % derivative of the muslce force constraints
    df_dlce = (dFce_dlce - dFse_dlce);
    df_dlce_opt = (dFce_dlce_opt - dFse_dlce_opt);
    df_dtheta0 = (dFce_dtheta0 - dFse_dtheta0);
    df_dW = (dFce_dW);
    df_ddlce = dFce_ddlce;
    df_dvmax_nor = dFce_dvmax_nor;
    df_dA = dFce_dA;
    df_dlce_slack_nor = dFce_dlce_slack_nor;
    df_dKp = dFce_dKp;
    % df_dDp = dFce_dDp;
    df_dgmax = dFce_dgmax;
    df_da = dFce_da;
    
    df_dKs = -dFse_dKs;
%     df_dlmt = -dFse_dlmt;
    df_dlt_slack = -dFse_dlt_slack;
    
end 