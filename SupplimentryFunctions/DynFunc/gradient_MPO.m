function grad = gradient_MPO(x, lmt, torque_m, act_m, d, M, N, S, P,...
                             muscle_par0, W1, W2, W3, W4, W5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The callback for calculating the gradient
%
% By: Huawei Wang
% Date: August 2, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    grad = zeros(size(x));
    
    % muscle parameters
    lce_opt = x(end - P + 1:end - P + M);
    lt_slack = x(end - P + M + 1:end - P + 2*M);
    theta0 = x(end - P + 2*M + 1:end - P + 3*M);
    Fmax = x(end - P + 3*M + 1:end);
    
    p_st = sum(N)*M*S;  % end of the state variables of all trials
    
    W1 = W1/sum(N);
    W2 = W2/sum(N);
    W3 = W3/sum(N);
    W4 = W4/sum(N);

    % gradient of the fit objective term
    % fit of muscle activations and joint torques
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
            [Fse, dFse_dlce, dFse_dlce_opt, dFse_dlt_slack, dFse_dtheta0] = ...
                tendenForce_Groote_diff(lmt_mea, lce, lce_opt, lt_slack, theta0);
            
            tor_opt = Fse.*Fmax*d_mea';
            
            grad(sta_st + sta_st_n + 2*M + 1: sta_st + sta_st_n + 3*M) = ...
                grad(sta_st + sta_st_n + 2*M + 1: sta_st + sta_st_n + 3*M) + ...
                2*W1*(tor_opt - tor_mea).*Fmax.*d_mea.*dFse_dlce;
            
            grad(p_st + 1:p_st + M) = grad(p_st + 1:p_st + M) + ...
                2*W1*(tor_opt - tor_mea).*Fmax.*d_mea.*dFse_dlce_opt;
            
            grad(p_st + M + 1:p_st + 2*M) = grad(p_st + M + 1:p_st + 2*M) + ...
                2*W1*(tor_opt - tor_mea).*Fmax.*d_mea.*dFse_dlt_slack;
            
            grad(p_st + 2*M + 1:p_st + 3*M) = grad(p_st + 2*M + 1:p_st + 3*M) + ...
                2*W1*(tor_opt - tor_mea).*Fmax.*d_mea.*dFse_dtheta0;
            
            grad(p_st + 3*M + 1:p_st + 4*M) = grad(p_st + 3*M + 1:p_st + 4*M) + ...
                2*W1*(tor_opt - tor_mea).*Fse.*d_mea;
            
            
            % calculate gradients
            grad(sta_st + sta_st_n + 4*M + 1:sta_st + sta_st_n + 5*M) = ...
                grad(sta_st + sta_st_n + 4*M + 1:sta_st + sta_st_n + 5*M) ... 
                + 2*W2*(act_opt - act_mea);
            
            % activation smoothness gradients
            if n < N(t)

                act_opt_nt = x(sta_st + sta_st_n + M*S + 4*M + 1: ...
                            sta_st + sta_st_n + M*S + 5*M);

                grad(sta_st + sta_st_n + M*S + 4*M + 1:sta_st + sta_st_n + M*S + 5*M) = ...
                grad(sta_st + sta_st_n + M*S + 4*M + 1:sta_st + sta_st_n + M*S + 5*M) ...
                + 2*W3*(act_opt_nt - act_opt);
            
                grad(sta_st + sta_st_n + 4*M + 1:sta_st + sta_st_n + 5*M) = ...
                grad(sta_st + sta_st_n + 4*M + 1:sta_st + sta_st_n + 5*M) ...
                - 2*W3*(act_opt_nt - act_opt);
            
                lce_nt = x(sta_st + sta_st_n + M*S + 2*M + 1: sta_st + sta_st_n + M*S + 3*M);
                lmt_mea_nt = lmt(mea_st_n + 1, :);
                [Fse_nt, dFse_nt_dlce, dFse_nt_dlce_opt, dFse_nt_dlt_slack, dFse_nt_dtheta0] = ...
                    tendenForce_Groote_diff(lmt_mea_nt, lce_nt, lce_opt, lt_slack, theta0);
                
                dsum_fse_dFse = -2*W4*(Fse_nt - Fse);
                dsum_fse_dFse_nt = 2*W4*(Fse_nt - Fse);
                
                % gradient of the muscle force changes
                grad(sta_st + sta_st_n + 2*M + 1: sta_st + sta_st_n + 3*M) = ...
                grad(sta_st + sta_st_n + 2*M + 1: sta_st + sta_st_n + 3*M) + ... 
                        dsum_fse_dFse.*dFse_dlce;

                % gradient of the muscle force changes
                grad(sta_st + sta_st_n + M*S + 2*M + 1: sta_st + sta_st_n + M*S + 3*M) = ...
                grad(sta_st + sta_st_n + M*S + 2*M + 1: sta_st + sta_st_n + M*S + 3*M) + ...
                        dsum_fse_dFse_nt.*dFse_nt_dlce;
                    
                % Fse derivatives
                grad(p_st + 1:p_st + M) = grad(p_st + 1:p_st + M) + ...
                dsum_fse_dFse.*dFse_dlce_opt;
            
                grad(p_st + M + 1:p_st + 2*M) = grad(p_st + M + 1:p_st + 2*M) + ...
                    dsum_fse_dFse.*dFse_dlt_slack;
                
                grad(p_st + 2*M + 1:p_st + 3*M) = grad(p_st + 2*M + 1:p_st + 3*M) + ...
                dsum_fse_dFse.*dFse_dtheta0;

                % Fse_nt derivatives
                grad(p_st + 1:p_st + M) = grad(p_st + 1:p_st + M) + ...
                    dsum_fse_dFse_nt.*dFse_nt_dlce_opt;
            
                grad(p_st + M + 1:p_st + 2*M) = grad(p_st + M + 1:p_st + 2*M) + ...
                    dsum_fse_dFse_nt.*dFse_nt_dlt_slack;
                
                grad(p_st + 2*M + 1:p_st + 3*M) = grad(p_st + 2*M + 1:p_st + 3*M) + ...
                    dsum_fse_dFse_nt.*dFse_nt_dtheta0;
                
            end
        end
    end

    grad(end - P + 1:end) = grad(end - P + 1:end) + ...
        2*W5*(x(end - P + 1:end)./muscle_par0 - 1)./muscle_par0/P;
    
end