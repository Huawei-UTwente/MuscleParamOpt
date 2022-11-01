function obj =  objective_MPO(x, lmt, torque_m, act_m, d, M, N, S, P,...
                          muscle_par0, W1, W2, W3, W4, W5, W6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The callback for calculating the objective
% Input: 
%        x: optimizing states [lce_opt, lt_slack, theta0]
%        num_nodes: direct collocation nodes (data nodes) in trajectory
%                   optimization
%        d:
%
% By: Huawei Wang
% Date: August 2, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    lce_opt = x(end - P + 1:end - P + M);
    lt_slack = x(end - P + M + 1:end - P + 2*M);
    theta0 = x(end - P + 2*M + 1:end - P + 3*M);
    Fmax = x(end - P + 3*M + 1:end);
    
    sum_tor = 0;
    sum_act = 0;
    sum_for = 0;  % very small term to eliminate  
    sum_act_smo = 0;
    sum_fse_smo = 0;
    
    for t = 1:length(N)
        sta_st = sum(N(1:t-1))*M*S;
        mea_st = sum(N(1:t-1));
    
        for n = 1:N(t)
            sta_st_n = (n-1)*M*S;
            mea_st_n = mea_st + n;
            
            % extract optimized states and measured data
            act_opt = x(sta_st + sta_st_n + 4*M + 1:sta_st + sta_st_n + 5*M);
            lce = x(sta_st + sta_st_n + 2*M + 1: sta_st + sta_st_n + 3*M);
            lmt_mea = lmt(mea_st_n, :);
            act_mea = act_m(mea_st_n, :);
            d_mea = d(mea_st_n, :);
            tor_mea = torque_m(mea_st_n, :);
            
            % calculate muscle force and joint torques
            Fse = tendenForce_Groote_MPO(lmt_mea, lce, lce_opt, lt_slack, theta0);
            tor_opt = Fse.*Fmax*d_mea';  
            
            % calculate the sum of square fit
            sum_tor = sum_tor + sum((tor_opt - tor_mea).^2);  % torque fit
            sum_act = sum_act + sum((act_opt - act_mea).^2);  % emg fit
            
            delta = 1e-6;
            sum_for = sum_for + sum(sqrt(Fse.^2 + delta));  % L1 regularization
            
            if n < N(t)
               act_opt_nt = x(sta_st + sta_st_n + M*S + 4*M + 1: ...
                              sta_st + sta_st_n + M*S + 5*M);
                          
               sum_act_smo = sum_act_smo + sum((act_opt_nt - act_opt).^2);
               
               lce_nt = x(sta_st + sta_st_n + M*S + 2*M + 1: sta_st + sta_st_n + M*S + 3*M);
               lmt_mea_nt = lmt(mea_st_n + 1, :);
               
               Fse_nt = tendenForce_Groote_MPO(lmt_mea_nt, lce_nt, lce_opt, lt_slack, theta0);
               sum_fse_smo = sum_fse_smo + sum((Fse_nt - Fse).^2);
               
            end
        end
    end
    
    obj_fit = (W1*sum_tor + W2*sum_act + W6*sum_for + W3*sum_act_smo + W4*sum_fse_smo)/sum(N);

    obj_par = sum((x(end - P + 1:end)./muscle_par0 - 1).^2)/P;
    
    obj = obj_fit + W5*obj_par;
    
end