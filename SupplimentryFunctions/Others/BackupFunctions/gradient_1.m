function grad = gradient_1(x, lmt, moment, u_meas, T, L, J, N, M, S, P, d,...
    W1, W2, W3, moment_range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The callback for calculating the gradient
%
% By: Huawei Wang
% Date: August 2, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    grad = zeros(size(x));

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
    index_u_matrix_max = find(u_matrix_max == u_matrix_re);

    if length(index_u_matrix_max) ~= M*T
        error('number of u_max index does not equal to the size of u_max')
    end
    
    u_meas_max = max(u_meas_re);

    for n = 1:N  % run over nodes        
        for t = 1:T  % run over trails
            for l = 1:L  % run over legs
                
                lce_ntl = x((n-1)*T*L*M*S + 2*T*L*M + (t-1)*L*M + (l-1)*M + 1: ...
                        (n-1)*T*L*M*S  + 2*T*L*M + (t-1)*L*M + l*M);  % all lce at the data node n
        
                u_ntl = x(N*T*L*M*S + (n-1)*T*L*M + (t-1)*L*M + (l-1)*M + 1: ...
                    N*T*L*M*S + (n-1)*T*L*M + (t-1)*L*M + l*M);  % all the muscle activations at the data node n
                
                lmt_ntl = lmt(n, (t-1)*L*M + (l-1)*M + 1: ...
                                 (t-1)*L*M + l*M);  % muscle tendon lengths of all muscles at the data node n 
                
                moment_ntl = moment(n, (t-1)*L*J + (l-1)*J + 1: ...
                                       (t-1)*L*J + l*J);  % experimental joint moments 
                             
                u_meas_ntl = u_meas(n, (t-1)*L*M + (l-1)*M + 1: ...
                                       (t-1)*L*M + l*M);  % experimental muscle activations
                                   
                moment_range_tl = moment_range((t-1)*L*J + (l-1)*J + 1:...
                    (t-1)*L*J + l*J);
                
                
                u_matrix_max_ntl = u_matrix_max((t-1)*M + 1:t*M);
                
                index_u_matrix_max_ntl = index_u_matrix_max((t-1)*M + 1:t*M);
                         
                index_u_matrix_max_ntl_recal = zeros(size(index_u_matrix_max_ntl));
                
                for m = 1:M
                    
                    ind_rl = index_u_matrix_max_ntl(m) - ((t-1)*M + m - 1)*N*L;
                    
                    if ind_rl == N
                        index_u_matrix_max_ntl_recal(m) = (N - 1)*L*T*M ...
                            + (t-1)*M*L + m;
                    elseif ind_rl == 2*N
                        index_u_matrix_max_ntl_recal(m) = (N - 1)*L*T*M ...
                            + (t-1)*M*L + M + m;
                    else
                        index_u_matrix_max_ntl_recal(m) = (rem(ind_rl, N) -...
                            1)*L*T*M + (t-1)*M*L + floor(ind_rl/N)*M + m;
                    end
                    
                end
                
                u_meas_max_ntl = u_meas_max((t-1)*M + 1:t*M);
                
                % calculate muscle force with the SEE tendon model
                [fse, dFse_dKs, dFse_dlce, dFse_dlce_opt, dFse_dtheta0,...
                    dFse_dlt_slack, dFse_dFmax] = ...
                tenden_force(lmt_ntl, lce_ntl, lce_opt, theta0, lt_slack, Ks, Fmax);
                
                % calculate joint moments based on the muscle force and moment arm
                mom = -fse*d';
                
                % calculate objective function, both moment and activation
                % are normalized by the experimental data.
                grad((n-1)*T*L*M*S + 2*T*L*M + (t-1)*L*M + (l-1)*M + 1: ...
                        (n-1)*T*L*M*S + 2*T*L*M + (t-1)*L*M + l*M) = ...
                        -2*W1*((moment_ntl - mom)./moment_range_tl.^2)*(-d).*dFse_dlce/(N*T);
                
                % generate grad vector
                grad(end - P + 2 + 1:end - P + 2 + M) = ...
                grad(end - P + 2 + 1:end - P + 2 + M) + ...
                        -2*W1*((moment_ntl - mom)./moment_range_tl.^2)*(-d).*dFse_dlce_opt/(N*T);

                grad(end - P + 2 + 10*M + 1:end) = ...
                grad(end - P + 2 + 10*M + 1:end) + ...
                        -2*W1*((moment_ntl - mom)./moment_range_tl.^2)*(-d).*dFse_dtheta0/(N*T);

                grad(end - P + 2 + 9*M + 1:end - P + 2 + 10*M) = ...
                grad(end - P + 2 + 9*M + 1:end - P + 2 + 10*M) + ...
                        -2*W1*((moment_ntl - mom)./moment_range_tl.^2)*(-d).*dFse_dlt_slack/(N*T);

                grad(end - P + 2 + 8*M + 1:end - P + 2 + 9*M) = ...
                grad(end - P + 2 + 8*M + 1:end - P + 2 + 9*M) + ...
                        -2*W1*((moment_ntl - mom)./moment_range_tl.^2)*(-d).*dFse_dKs/(N*T);
                
                grad(end - P + 2 + 5*M + 1:end - P + 2 + 6*M) = ...
                grad(end - P + 2 + 5*M + 1:end - P + 2 + 6*M) + ...
                        -2*W1*((moment_ntl - mom)./moment_range_tl.^2)*(-d).*dFse_dFmax/(N*T);
                
                grad(N*T*L*M*S + (n-1)*T*L*M + (t-1)*L*M + (l-1)*M + 1: ...
                    N*T*L*M*S + (n-1)*T*L*M + (t-1)*L*M + l*M) = ...
                grad(N*T*L*M*S + (n-1)*T*L*M + (t-1)*L*M + (l-1)*M + 1: ...
                    N*T*L*M*S + (n-1)*T*L*M + (t-1)*L*M + l*M) ...
                        -2*W2*(u_meas_ntl./u_meas_max_ntl -...
                        u_ntl./u_matrix_max_ntl)./u_matrix_max_ntl/(N*T);
                
                grad(N*T*L*M*S + index_u_matrix_max_ntl_recal) = ...
                grad(N*T*L*M*S + index_u_matrix_max_ntl_recal) ...
                        +2*W2*(u_meas_ntl./u_meas_max_ntl -...
                        u_ntl./u_matrix_max_ntl).*u_ntl./u_matrix_max_ntl.^2/(N*T);
                
                % objective function of the absolute muscle activation change.
                if n > 1
                    
                    u_ntl_1 = x(N*T*L*M*S + (n-2)*T*L*M + (t-1)*L*M + (l-1)*M + 1: ...
                    N*T*L*M*S + (n-2)*T*L*M + (t-1)*L*M + l*M);  % all the muscle activations at the data node n-1
                
                    grad(N*T*L*M*S + (n-1)*T*L*M + (t-1)*L*M + (l-1)*M + 1: ...
                    N*T*L*M*S + (n-1)*T*L*M + (t-1)*L*M + l*M) = ...
                    grad(N*T*L*M*S + (n-1)*T*L*M + (t-1)*L*M + (l-1)*M + 1: ...
                    N*T*L*M*S + (n-1)*T*L*M + (t-1)*L*M + l*M) ...
                        + 2*W3*(u_ntl - u_ntl_1)/(N*T);

                    grad(N*T*L*M*S + (n-2)*T*L*M + (t-1)*L*M + (l-1)*M + 1: ...
                    N*T*L*M*S + (n-2)*T*L*M + (t-1)*L*M + l*M) = ...
                    grad(N*T*L*M*S + (n-2)*T*L*M + (t-1)*L*M + (l-1)*M + 1: ...
                    N*T*L*M*S + (n-2)*T*L*M + (t-1)*L*M + l*M) ...
                        -2*W3*(u_ntl - u_ntl_1)/(N*T);
                end
                
            end
        end
    end
end