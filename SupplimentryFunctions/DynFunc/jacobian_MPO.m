function Jac = jacobian_MPO(x, M, S, C, N, P, lmt, t_em, hs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The callback for calculating the Jacobian
%
% By: Huawei Wang
% Date: August 3, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % muscle parameters
    par = x(end - P + 1:end);
 
    % Initilize the empty nonzero Jacobian vector
    Jac = [];

    for t = 1:length(N)      % Run through all data trials
        
        hs_t = hs(t);        % Time interval between the data nodes
        n_em = t_em(t)/hs_t; % Number of data nodes that caused by electromechanical delay (_em)
        
        x_st = sum(N(1:t-1))*M*S; % The index of the first state parameter of current data trial 
        lmt_st = sum(N(1:t-1));   % The index of the first lmt data of current data trial
        
        % Looking back to the data nodes that caused by _em and _rf delays 
        stn_em1 = ceil(n_em)*M*S;
        stn_em2 = floor(n_em)*M*S;
        
        for n = 1: N(t)-1   % Run through all data nodes (using middle point method)
            
            x_stn = x_st + (n-1)*M*S; 	% The index right before the first state of current data node
            x_stn1 = x_st + n*M*S;      % The index of the last state of current data node
            x_stn2 = x_st + (n+1)*M*S;  % The index of last state of next data node
            
            % Calculate the x, dx, and lmt based on the middle point method
            x_tn = (x(x_stn1+1: x_stn2) + x(x_stn+1 : x_stn1))/2;
            dx_tn = (x(x_stn1+1: x_stn2) - x(x_stn+1 : x_stn1))./hs_t;
            
            lmt_tn = (lmt(lmt_st + n, :) + lmt(lmt_st + n + 1, :))/2; 
            
            % If the number of current data node is smaller than the data node
            % required by the _em delay. Then the activation
            % dynamic and reflex control constraints are not included,
            % since there are no state can be used to formulate these
            % feedbacks.
            if n <= ceil(n_em)
                
                % directCollocationDyn_diff_RPO1 is the derivative function
                % for directCollocationDyn_RPO1. Only df_dx, and fx_ddx are
                % calculate, since x_em and x_rf are not available
                [df_dx_tn, df_ddx_tn, df_dp] = ...
                directCollocationDyn_diff_MPO1(x_tn, dx_tn, ...
                    lmt_tn, par, M, P, S);
                
                % (middle point method) assign back derivatives to data nodes
                df_dx1 = df_dx_tn./2 - df_ddx_tn./hs_t;
                df_dx2 = df_dx_tn./2 + df_ddx_tn./hs_t;

                % extract nonzero elements in the jacobian matrix, row by
                % row
                for r = 1:M*(C-1)
                    Jac = [Jac, nonzeros([df_dx1(r, :), df_dx2(r, :), df_dp(r, :)])'];
                end
                
            % Else if the number of current data node is larger than the
            % _em, Then the activation dynamics is included      
            else
            
                % If n_em is not integer, the _em delay data nodes are
                % calculated using linear interpolation between the near by
                % data nodes. If n_em is integer, then directly find the
                % delayed data nodes.
                if n_em == floor(n_em)
                    w1_em = 0;
                    w2_em = 1;
                else
                    w1_em = n_em - floor(n_em);
                    w2_em = ceil(n_em) - n_em;
                end
                
                % _em delayed data node of x1
                x_em1 = x(x_stn-stn_em1+1 : x_stn1-stn_em1).*w1_em  ...
                        + x(x_stn-stn_em2+1 : x_stn1-stn_em2).*w2_em;
                
                % _em delayed data node of x2
                x_em2 = x(x_stn1-stn_em1+1 : x_stn2-stn_em1).*w1_em  ...
                        + x(x_stn1-stn_em2+1 : x_stn2-stn_em2).*w2_em;
                
                % Middle point method to get the _em delayed variable for
                % the activation dynamics
                x_em_tn = (x_em1 + x_em2)/2;
               
                % directCollocationDyn_diff_RPO2 is the derivative function
                % for directCollocationDyn_RPO2. df_dx_em is included here,
                % however, df_dx_rf is still not included, since x_rf is
                % not available yet.
                [df_dx_tn, df_ddx_tn, df_dx_em_tn, df_dp] = ...
                directCollocationDyn_diff_MPO2(x_tn, dx_tn, x_em_tn, ...
                    lmt_tn, par, M, P, S);
                
                % (middle point method) assign back derivatives to data nodes
                df_dx1 = df_dx_tn./2 - df_ddx_tn./hs_t;
                df_dx2 = df_dx_tn./2 + df_ddx_tn./hs_t;

                df_dx_em1 = df_dx_em_tn.*w1_em/2;
                df_dx_em2 = df_dx_em_tn.*(w1_em + w2_em)/2;
                df_dx_em3 = df_dx_em_tn.*w2_em/2;

                % extract nonzero elements in the jacobian matrix, row by
                % row
                for r = 1:M*C
                    if r <= 2*M
                        
                        Jac = [Jac, nonzeros([df_dx1(r, :), df_dx2(r, :)])'];

                    elseif r <= 3*M
                        
                        df_dx_x_em = jacDelay_em(n_em, df_dx_em1(2*M+1:3*M, :),...
                            df_dx_em2(2*M+1:3*M, :), df_dx_em3(2*M+1:3*M, :),...
                            df_dx1(2*M+1:3*M, :), df_dx2(2*M+1:3*M, :), M, S);
                         
                        Jac = [Jac, nonzeros(df_dx_x_em(r - 2*M, :))'];
  
                    else
                        
                        Jac = [Jac, nonzeros([df_dx1(r, :), df_dx2(r, :), df_dp(r, :)])'];

                    end
                end
            end
        end
    end
end