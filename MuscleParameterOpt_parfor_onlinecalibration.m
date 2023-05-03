%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is to optimize muscle parameters from experimental data.
%
% The experimental data should contain joint angles, joint torques, and
% muscle activations (EMG).
%
% Four muscle parameters (l_opt, lt_slack, theta0, and Fmax) will be
% optimized for each muscle.
%
% Current code is for the ankle joint muscles including: 
% M/L Gastrocneius, Soleus; Tibiais Anterior. 
%
% However, it is easy to extend to other muscles and joints. Just change
% the inputs in this script, no need to go deeper functions. The deeper
% functions were coded in a general way that can handle any number of
% muscles and joints.
%
% Gradient based optimizer is used here to find the best muscle parameters that
% can explain the experimental data.
%
% By: Huawei Wang
% Date: Augst 1, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

homeDataPath = 'exampleData';
trialNamesOpt = ["lift_0"];
muscleNames = ["soleus_r","gas_r"];
coordinateNames = ["ankle_angle_r"];

MuscleParameterOpt_parfor_func(homeDataPath, trialNamesOpt, ...
                                        muscleNames, coordinateNames)

function MuscleParameterOpt_parfor_func(homeDataPath, trialNamesOpt, ...
                                        muscleNames, coordinateNames)
               
    T = length(trialNamesOpt);  % number of data trails involved in the optimization, each data trial
            % could have different data nodes.
    
    t_em = 0.1 + zeros(1, T);  % electrialmechaincal delay inside muscle physiologies, around
                               % 100ms. Could be different among data trials
                               
    J = length(coordinateNames);  % number of joint in each side of leg
    M = length(muscleNames);  % number of muscle in each side of leg
    S = 5;  % number of states of each muscle model
    C = 4;  % number of constraints of each muscle model

    hs = zeros(1, T);

    for t = 1:T  % load trial data and get muscle parameters from the averaged gait data
        
        trial = trialNamesOpt(t);

        % get muscle parameters
        

        % get muscle activation
        act_dataset = importdata(sprintf('%s/%s/Activations.sto', homeDataPath, trial));
        act = act_dataset.data(100:600, 2:3);
        
        % get joint torques
        tor_dataset = importdata(sprintf('%s/%s/id.sto', homeDataPath, trial));
        tor = tor_dataset.data(100:600, 2);
        
        % get muscle lengths
        FiberLen_dataset = importdata(sprintf('%s/%s/FibreLengths.sto', homeDataPath, trial));
        TendonLen_dataset = importdata(sprintf('%s/%s/TendonLength.sto', homeDataPath, trial));
        muslen = FiberLen_dataset.data(100:600, 2:3) + TendonLen_dataset.data(100:600, 2:3);

        % get moment arms
        ma_dataset = importdata(sprintf('%s/%s/ma_ankle_angle_r.sto', homeDataPath, trial));
        ma = ma_dataset.data(100:600, 2:3);
        
        % get the time interval
        hs(t) = 0.1;
        
        homeSavingPath = sprintf('%s/%s/opt', homeDataPath, trial);

    end
    
    N = length(torque);

    % initial muscle parameters from the initial xml
    mus_par0 = importdata(sprintf('%s/mus_param.txt', homeDataPath)); 
    P = length(mus_par0);
    
    % optimization saving path
    mkdir(homeSavingPath);
    
    %% define default muscle parameters
    % weights of the terms of the objective function
    W1 = 1;  % weight of joint torque fits
    W2 = 2;  % weight of muscle activation fits
    W3 = 0.1;  % weight of muscle activation smoothness
    W4 = 0.1;  % weight of muscle force smoothness
    W5 = 0.01; % weight of diversity of the optimizing parameters
    W6 = 0.005;  % very small number to eliminate co-contraction
    
    Prange = 0.5;
    
    %% optimization boundaries and setups
    x_lb1 = [zeros(sum(N), M) + 0.001,...              % activation
            zeros(sum(N), M) - 30,...                  % d_activation
            zeros(sum(N), M) + mus_par0(1:M)*0.35,...    % lce
            zeros(sum(N), M) - mus_par0(1:M)*15, ....   % dlce
            zeros(sum(N), M) + 0.001];                 % nerual stimulation

    x_lb2 = reshape(x_lb1', [1, sum(N)*M*S]);
    
    % l_opt, lt_slack have range of Prange, penn_ang is locked, Fmax has
    % range of 2*Prange
    par_rf_lb = [mus_par0(1:2*M)*(1 - Prange), mus_par0(2*M + 1:3*M), ...
                 mus_par0(3*M + 1:4*M)*(1 - Prange)];

    x_lb3 = [x_lb2, par_rf_lb];

    % set up optimizing parameter upper bounds.  
    % lce can be maximum 2 time the lce_opt
    x_ub1 = [zeros(sum(N), M) + 1,...                  % activation
             zeros(sum(N), M) + 30,...                 % d_activation
             zeros(sum(N), M) + mus_par0(1:M)*1.8,...   % lce
             zeros(sum(N), M) + mus_par0(1:M)*15, ...   % dlce
             zeros(sum(N), M) + 1];                    % nerual stimulation

    x_ub2 = reshape(x_ub1', [1, sum(N)*M*S]);

    % l_opt, lt_slack have range of Prange, penn_ang is locked, Fmax has
    % range of 2*Prange
    par_rf_ub = [mus_par0(1:2*M)*(1 + Prange), mus_par0(2*M + 1:3*M), ...
                 mus_par0(3*M + 1:4*M)*(1 + Prange*3)]; 

    x_ub3 = [x_ub2, par_rf_ub];
    
    %% do optimizations with ipopt
    
    % Set the IPOPT options.
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy           = 'adaptive';
    options.ipopt.tol                   = 1e-3;
    options.ipopt.max_iter              = 5000;
    options.ipopt.linear_solver         = 'mumps';
    % options.ipopt.max_cpu_time         = 1000;
    
    % optimization parameters
    auxdata.mus_par0 = mus_par0;
    auxdata.t_em = t_em;
    
    auxdata.M = M;
    auxdata.S = S;
    auxdata.C = C;
    auxdata.N = N;
    auxdata.T = T;
    auxdata.J = J;
    auxdata.P = P;
    auxdata.lmt = lmt;
    auxdata.hs = hs;
    auxdata.torque = torque;
    auxdata.mus_act = mus_act;
    auxdata.d = d;
    
    % weight of each objective function term
    auxdata.W1 = W1;
    auxdata.W2 = W2;
    auxdata.W3 = W3;
    auxdata.W4 = W4;
    auxdata.W5 = W5;
    auxdata.W6 = W6;
    
    % jacobian structure
    [row, col] = jacobianstructure_ipopt_rc_MPO(auxdata);
    auxdata.row = row;
    auxdata.col = col;

    % boundaries for the state parameters
    options.lb = x_lb3;
    options.ub = x_ub3;
    
    rCons = M*C*sum(N - ceil(t_em./hs) - 1) ...
                + M*(C-1)*sum(ceil(t_em./hs));
            
    auxdata.rJac = rCons;
    auxdata.cJac = length(x_ub3);
            
    % boundaries for constraints, all zeros means all of them are equality
    % constraints.
    options.cl = zeros(1, rCons);
    options.cu = zeros(1, rCons);
    
    % Set the IPOPT options.
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy           = 'adaptive';
    options.ipopt.tol                   = 1e-3;
    options.ipopt.max_iter              = 5000;
    options.ipopt.linear_solver         = 'mumps';

    % The callback functions.
    funcs.objective         = @(x) objective_ipopt_MPO(x, auxdata);
    funcs.constraints       = @(x) constraints_ipopt_MPO(x, auxdata);
    funcs.gradient          = @(x) gradient_ipopt_MPO(x, auxdata);
    funcs.jacobian          = @(x) jacobian_ipopt_MPO(x, auxdata);
    funcs.jacobianstructure = @() jacobianstructure_ipopt_MPO(auxdata);

    auxdata.options = options;
    auxdata.funcs = funcs;

    auxdata.folder = homeSavingPath;
    
    % run optimizations in parallel 
    parfor opt = 1:10
        do_optimization_MPO(opt, auxdata)
    end
    
end
    