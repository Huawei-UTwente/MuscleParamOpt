function [row, col] = jacobianStructure_MPO(M, S, C, N, P, lmt, t_em, hs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sepcify the structure of Jacobian,
% in case the size of Jacobian is too large. 
%
% By: Huawei Wang
% Date: August 3, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% number of repetitions in generate the jac matrix. Normally larger
    % than 1, in case some location get false zero due to the random
    % selected states x.
    num_test = 2;  

    % Initilize the empty row and col vector for nonzero Jacobian vector
    row = [];
    col = [];

    % Get the controller parameters, reshape the fse, lce reflex control 
    % gains from vector to matrix
    
    p_st = sum(N)*M*S; % The index of the last state parameter of all data trial 
    
    for t = 1:length(N)       % Run through all data trials
        
        hs_t = hs(t);         % Time interval between the data nodes
        n_em = t_em(t)/hs_t;  % Number of data nodes that caused by electromechanical delay (_em)
        
        x_st = sum(N(1:t-1))*M*S; % The index of the first state parameter of current data trial 
        lmt_st = sum(N(1:t-1));   % The index of the first lmt data of current data trial
        
        % Looking back to the data nodes that caused by _em and _rf delays 
        stn_em1 = ceil(n_em)*M*S;
        % stn_em2 = floor(n_em)*M*S;

        % Index of the first constraint vector of the current data trial
        cons_st = M*C*sum(N(1:t-1) - ceil(t_em(1:t-1)./hs(1:t-1)) - 1) ...
            + M*(C-1)*sum(ceil(t_em(1:t-1)./hs(1:t-1)));
        
        % Number of constraint vectors that do not have either _em and _rf
        % be aware that they have different number of constraints in each 
        % secinarios. With neither _em and _rf, number of constraints is
        % C - 2; with only _em, number of constraints is C - 1.
        cons_st_em = ceil(n_em)*M*(C-1);
              
        for n = 1 : N(t)-1  % Run through all data nodes (using middle point method)
            
            x_stn = x_st + (n-1)*M*S; % The index right before the first state of current data node
            % x_stn1 = x_st + n*M*S;      % The index of the last state of current data node
            % x_stn2 = x_st + (n+1)*M*S;  % The index of last state of next data node
                     
            % If the number of current data node is smaller than the data node
            % required by the _em delay. Then the activation
            % dynamic and reflex control constraints are not included,
            % since there are no state can be used to formulate these
            % feedbacks.
            if n <= ceil(n_em)
                
                df_dx1 = zeros(M*(C-1), M*S);
                df_dx2 = zeros(M*(C-1), M*S);
                df_dp = zeros(M*(C-1), P);
                for test = 1:num_test   % run multiple repetitions
                    
                    % generate randomized x1 and x2, instead of extracting
                    % from the optimizing variables x. Randomization can
                    % largly avoid false zero elements.
                    x1 = rand(1, M*S);
                    x2 = rand(1, M*S);
                    
                    x1(2*M+1:3*M) = x1(2*M+1:3*M)*0.01+0.05;
                    x1(3*M+1:4*M) = -x1(3*M+1:4*M)*0.2+ 0.1;

                    x2(2*M+1:3*M) = x2(2*M+1:3*M)*0.01+0.06;
                    x2(3*M+1:4*M) = -x2(3*M+1:4*M)*0.3 + 0.2;
                    
                    par = rand(1, P);  % randomize muscle parameters
            
                    par(1:M) = 0.04 + 0.02*par(1:M);  % lce_opt
                    par(M+1:2*M) = par(M+1:2*M)*0.1 + 0.3;  % lt_slack
                    par(2*M+1:3*M) = par(2*M+1:3*M)*pi/18 + pi/18;  % theta0
                    par(3*M+1:4*M) = par(3*M+1:4*M)*50 + 10;  % Fmax

                    % Calculate the x, dx, and lmt based on the middle point method
                    x_tn = (x1 + x2)/2;
                    dx_tn = (x1 - x2)./hs_t;

                    lmt_tn = (lmt(lmt_st + n, :) + lmt(lmt_st + n + 1, :))/2; 
                    
                    % directCollocationDyn_diff_RPO1 is the derivative function
                    % for directCollocationDyn_RPO1. Only df_dx, and fx_ddx are
                    % calculate, since x_em and x_rf are not available
                    [df_dx_tn, df_ddx_tn, df_dpar] = ...
                    directCollocationDyn_diff_MPO1(x_tn, dx_tn, ...
                        lmt_tn, par, M, P, S);

                    % (middle point method) assign back derivatives to data nodes
                    df_dx1 = df_dx1 + df_dx_tn./2 - df_ddx_tn./hs_t;
                    df_dx2 = df_dx2 + df_dx_tn./2 + df_ddx_tn./hs_t;
                    
                    df_dp = df_dp + df_dpar;
                    
                end
                    
                % extract the indexes of nonzero elements in the jacobian 
                % matrix, row by row
                row_st = cons_st + (n - 1)*M*(C-1);

                for r = 1:M*(C-1)
                    if r <= 2*M
                        [row_i, col_i] = find([df_dx1(r, :), df_dx2(r, :)]);

                        row_i = row_i + row_st + r - 1;
                        col_i = col_i + x_stn;

                        row = [row, row_i];
                        col = [col, col_i];
                    
                    else
                        [row_i, col_i] = find([df_dx1(r, :), df_dx2(r, :)]);

                        row_i = row_i + row_st + r - 1;
                        col_i = col_i + x_stn;
                        
                        [row_ip, col_ip] = find(df_dp(r, :));
                        
                        row_ip = row_ip + row_st + r - 1;
                        col_ip = col_ip + p_st;
                        
                        row = [row, row_i, row_ip];
                        col = [col, col_i, col_ip];
                    end
                end        

            % Else if the number of current data node is larger than the
            % _em, but smaller than _em + _rf. Then the activation dynamics
            % is included, the reflex control constraint will not be 
            % included.
            else
                
                df_dx1 = zeros(M*C, M*S);
                df_dx2 = zeros(M*C, M*S);

                df_dx_em1 = zeros(M*C, M*S);
                df_dx_em2 = zeros(M*C, M*S);
                df_dx_em3 = zeros(M*C, M*S);
                
                df_dp = zeros(M*C, P);
                
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
                   
                for test = 1:num_test
                    
                    % generate randomized x1, x2, and x_em, instead of extracting
                    % from the optimizing variables x. Randomization can
                    % largly avoid false zero elements.
                    x1 = rand(1, M*S);
                    x2 = rand(1, M*S);

                    x1(2*M+1:3*M) = x1(2*M+1:3*M)*0.01+0.05;
                    x1(3*M+1:4*M) = -x1(3*M+1:4*M)*0.2+ 0.1;

                    x2(2*M+1:3*M) = x2(2*M+1:3*M)*0.01+0.06;
                    x2(3*M+1:4*M) = -x2(3*M+1:4*M)*0.3 + 0.2;
                    
                    par = rand(1, P);  % randomize muscle parameters
            
                    par(1:M) = 0.04 + 0.02*par(1:M);  % lce_opt
                    par(M+1:2*M) = par(M+1:2*M)*0.1 + 0.3;  % lt_slack
                    par(2*M+1:3*M) = par(2*M+1:3*M)*pi/18 + pi/18;  % theta0
                    par(3*M+1:4*M) = par(3*M+1:4*M)*50 + 10;  % Fmax

                    % calculate the x and dx based on the middle point method
                    x_tn = (x1 + x2)/2;
                    dx_tn = (x1 - x2)./hs_t;

                    lmt_tn = (lmt(lmt_st + n, :) + lmt(lmt_st + n + 1, :))/2; 
                    
                    x_em_tn = rand(1, M*S);

                    % directCollocationDyn_diff_RPO2 is the derivative function
                    % for directCollocationDyn_RPO2. df_dx_em is included here,
                    % however, df_dx_rf is still not included, since x_rf is
                    % not available yet.
                    [df_dx_tn, df_ddx_tn, df_dx_em_tn, df_dpar] = ...
                    directCollocationDyn_diff_MPO2(x_tn, dx_tn, x_em_tn, ...
                        lmt_tn, par, M, P, S);
                    
                    % (middle point method) assign back derivatives to data nodes
                    df_dx1 = df_dx1 + df_dx_tn./2 - df_ddx_tn./hs_t;
                    df_dx2 = df_dx2 + df_dx_tn./2 + df_ddx_tn./hs_t;

                    df_dx_em1 = df_dx_em1 + df_dx_em_tn.*w1_em/2;
                    df_dx_em2 = df_dx_em2 + df_dx_em_tn/2;
                    df_dx_em3 = df_dx_em3 + df_dx_em_tn*w2_em/2;
                    
                    df_dp = df_dp + df_dpar;

                end
                    
                % extract nonzero elements in the jacobian matrix, row by
                % row
                row_st = cons_st + cons_st_em + (n - ceil(n_em) - 1)*M*C;

                for r = 1:M*C
                    if r <= 2*M  % derivative constraints
                        
                        [row_i, col_i] = find([df_dx1(r, :), df_dx2(r, :)]);

                        row_i = row_i + row_st + r - 1;
                        col_i = col_i + x_stn;
                        
                        row = [row, row_i];
                        col = [col, col_i];

                    elseif r <= 3*M  % activation constraints
                        
                        df_dx_x_em = jacDelay_em(n_em, df_dx_em1(2*M+1:3*M, :),...
                            df_dx_em2(2*M+1:3*M, :), df_dx_em3(2*M+1:3*M, :),...
                            df_dx1(2*M+1:3*M, :), df_dx2(2*M+1:3*M, :), M, S);
                         
                        [row_i, col_i] = find(df_dx_x_em(r - 2*M, :));
                                
                        row_i = row_i + row_st + r - 1;
                        col_i = col_i + x_stn - stn_em1;
                        
                        row = [row, row_i];
                        col = [col, col_i];
                        
                    elseif r <= 4*M  % muscle dynamic constraints
                        
                        [row_i, col_i] = find([df_dx1(r, :), df_dx2(r, :)]);

                        row_i = row_i + row_st + r - 1;
                        col_i = col_i + x_stn;
                        
                        [row_ip, col_ip] = find(df_dp(r, :));

                        row_ip = row_ip + row_st + r - 1;
                        col_ip = col_ip + p_st;
                        
                        row = [row, row_i, row_ip];
                        col = [col, col_i, col_ip];
                        
                    end 
                end
            end
        end
    end
end