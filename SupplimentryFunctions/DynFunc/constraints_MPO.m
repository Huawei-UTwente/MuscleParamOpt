function cons = constraints_MPO(x, M, S, C, N, P, lmt, t_em, hs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The callback for calculating the constraints
%
% By: Huawei Wang
% Date: August 3, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initilize the constraint vector size. Due to the _em and _rf delay,
    % not every data node has the constraint, only exist when the closed-loop
    % system with time delay can be fomulated.
    
    % Number of dynamic constraints of the optimization problem, because of
    % the reflex and electromechanical delay, the totoal constraints are
    % not M*C*sum(N), but smaller:
    rCons = M*C*sum(N - 1 - ceil(t_em./hs)) ...  % number of full elements constraints
            + M*(C-1)*sum(ceil(t_em./hs)); % number of constrinats that does not have activation dynamics
        
    % Initilize the constraint vector
    cons = zeros(1, rCons);
    
    % Extract the muscle parameters parameters 
    par = x(end - P + 1:end);

    for t = 1:length(N)      % Run through all data trials
        
        hs_t = hs(t);        % Time interval between the data nodes
        n_em = t_em(t)/hs_t; % Number of data nodes that caused by electromechanical delay (_em)
        
        x_st = sum(N(1:t-1))*M*S; % The index of the first state parameter of current data trial 
        lmt_st = sum(N(1:t-1));   % The index of the first lmt data of current data trial
        
        % Looking back to the data nodes that caused by _em and _rf delays 
        stn_em1 = ceil(n_em)*M*S;
        stn_em2 = floor(n_em)*M*S;

        % Index of the first constraint vector of the current data trial
        cons_st = M*C*sum(N(1:t-1) - ceil(t_em(1:t-1)./hs(1:t-1)) - 1) ...
            + M*(C-1)*sum(ceil(t_em(1:t-1)./hs(1:t-1)));
        
        % Number of constraint vectors that do not have either _em and _rf
        % be aware that they have different number of constraints in each 
        % secinarios. With neither _em and _rf, number of constraints is
        % C - 2; with only _em, number of constraints is C - 1.
        cons_st_em = ceil(n_em)*M*(C-1);
        
        for n = 1:N(t)-1    % Run through all data nodes (using middle point method)
            
            x_stn = x_st + (n-1)*M*S;   % The index right before the first state of current data node
            x_stn1 = x_st + n*M*S;      % The index of the last state of current data node
            x_stn2 = x_st + (n+1)*M*S;  % The index of last state of next data node
            
            % Calculate the x, dx, and lmt based on the middle point method
            x_tn = (x(x_stn1 + 1 : x_stn2) + x(x_stn + 1 : x_stn1))/2;
            dx_tn = (x(x_stn1 + 1 : x_stn2) - x(x_stn + 1 : x_stn1))./hs_t;
            
            lmt_tn = (lmt(lmt_st + n, :) + lmt(lmt_st + n + 1, :))/2; 
            
            % If the number of current data node is smaller than the data node
            % required by the _em delay. Then the activation
            % dynamic and reflex control constraints are not included,
            % since there are no state can be used to formulate these
            % feedbacks.
            if n <= ceil(n_em)
                
               % directCollocationDyn_RPO1 is the function that does not
               % constain activation dynamic and reflex control constraints
               cons(cons_st + (n-1)*M*(C-1)+1:cons_st + n*M*(C-1)) =...
                   directCollocationDyn_MPO1(x_tn, dx_tn, lmt_tn, par, M);
               
            % Else if the number of current data node is larger than the
            % _em, but smaller than _em + _rf. Then the activation dynamics
            % is included, the reflex control constraint will not be 
            % included.
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
                    
                % directCollocationDyn_RPO2 is the function that contains
                % activation dynamics, but does not contain reflex control
                % constraints
                cons(cons_st + cons_st_em + (n-ceil(n_em)-1)*M*C+1 ...
                    :cons_st + cons_st_em + (n-ceil(n_em))*M*C) =...
                   directCollocationDyn_MPO2(x_tn, dx_tn, x_em_tn, lmt_tn,...
                                             par, M);
                                         
            end
        end
    end
end