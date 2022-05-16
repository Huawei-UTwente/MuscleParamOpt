function obj =  objective(x, lmt, moment, u_meas, T, L, J, N, M, S, P, d,...
                          W1, W2, W3, moment_range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The callback for calculating the objective
% Input: 
%        x: optimizing states [alpha, beta1, beta2, Am, Lce0, vmax, gmax,
%                              W, A, Fmax, Kp, Dp, PEE_slack, Ks, SEE_slack,
%                              theta0]
%        num_nodes: direct collocation nodes (data nodes) in trajectory
%                   optimization
%        d:
%
% By: Huawei Wang
% Date: August 2, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    lce_opt = x(end - P + 2 + 1:end - P + 2 + M);
    theta0 = x(end - P + 2 + 10*M + 1:end);
    lt_slack = x(end - P + 2 + 9*M + 1:end - P + 2 + 10*M);
    Ks = x(end - P + 2 + 8*M + 1:end - P + 2 + 9*M);
    Fmax = x(end - P + 2 + 5*M + 1:end - P + 2 + 6*M);
    
    u_matrix = reshape(x(N*T*M*L*S + 1:N*T*M*L*S + N*T*M*L), T*M*L, N)';
    u_matrix_re = zeros(N*L, T*M);
    u_meas_re = zeros(N*L, T*M);
    
    if L == 1
        u_matrix_re = u_matrix;
        u_meas_re = u_meas;
    elseif L == 2
        for t = 1:T
            u_matrix_re(:, (t-1)*M + 1:t*M) = [u_matrix(:, (t-1)*M*L + 1:(t-1)*M*L + M);...
                u_matrix(:, (t-1)*M*L + M + 1:(t-1)*M*L + 2*M)];
            
            u_meas_re(:, (t-1)*M + 1:t*M) = [u_meas(:, (t-1)*M*L + 1:(t-1)*M*L + M);...
                u_meas(:, (t-1)*M*L + M + 1:(t-1)*M*L + 2*M)];
            
        end
    end
    
    u_matrix_max = max(u_matrix_re);
   
    u_meas_max = max(u_meas_re);

    obj_fm = 0;
    obj_fu = 0;
    obj_fdu = 0;
    
    for n = 1:N  % run over nodes        
        for t = 1:T  % run over trails
            for l = 1:L  % run over legs
                
                lce_ntl = x((n-1)*T*L*M*S + 2*T*L*M + (t-1)*L*M + (l-1)*M + 1: ...
                            (n-1)*T*L*M*S + 2*T*L*M + (t-1)*L*M + l*M);  % all lce at the data node n
        
                u_ntl = x(N*T*L*M*S + (n-1)*T*L*M + (t-1)*L*M + (l-1)*M + 1: ...
                    N*T*L*M*S + (n-1)*T*L*M + (t-1)*L*M + l*M);  % all the muscle activations at the data node n
                
                lmt_ntl = lmt(n, (t-1)*L*M + (l-1)*M + 1: ...
                                 (t-1)*L*M + l*M);  % muscle tendon lengths of all muscles at the data node n 
                
                moment_ntl = moment(n, (t-1)*L*J + (l-1)*J + 1: ...
                                       (t-1)*L*J + l*J);
                             
                u_meas_ntl = u_meas(n, (t-1)*L*M + (l-1)*M + 1: ...
                                       (t-1)*L*M + l*M);
                                   
                moment_range_tl = moment_range((t-1)*L*J + (l-1)*J + 1:...
                    (t-1)*L*J + l*J);
                
                
                u_matrix_max_ntl = u_matrix_max((t-1)*M + 1:t*M);
                                            
                u_meas_max_ntl = u_meas_max((t-1)*M + 1:t*M);
                                   
                % calculate muscle force with the SEE tendon model
                [fse, ~, ~, ~, ~, ~] = ...
                tenden_force(lmt_ntl, lce_ntl, lce_opt, theta0, lt_slack,...
                             Ks, Fmax);
                     
                % calculate joint moments based on the muscle force and
                % moment arm
                mom = -fse*d';
                
                % calculate objective function, both moment and activation
                % are normalized by the experimental data.
                obj_fm = obj_fm + sum(((moment_ntl - mom)./moment_range_tl).^2)/(N*T);
                
                obj_fu = obj_fu + sum((u_meas_ntl./u_meas_max_ntl-...
                    u_ntl./u_matrix_max_ntl).^2)/(N*T);
                
                % objective function of the absolute muscle activation
                % change.
                if n > 1
                    u_ntl_1 = x(N*T*L*M*S + (n-2)*T*L*M + (t-1)*L*M + (l-1)*M + 1: ...
                    N*T*L*M*S + (n-2)*T*L*M + (t-1)*L*M + l*M);  % all the muscle activations at the data node n-1

                    obj_fdu = obj_fdu + sum((u_ntl - u_ntl_1).^2)/(N*T);
                end 
            end
        end
    end
    
    obj = (W1*obj_fm + W2*obj_fu + W3*obj_fdu);
    
end