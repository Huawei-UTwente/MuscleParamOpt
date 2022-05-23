%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative test of the contaction dynamics
%
% By: Huawei Wang
% Date: 12/06/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

subj = 6;

% specify the data trials
trialNames = ["walk_09", "walk_18", "walk_27", "walk_45", "walk_54"];

% initlize parameters
T = 3;                      % number of data trials
N = 50 + zeros(1, T);       % number of data node at each data trial
t_em = 0.1 + zeros(1, T);   % muscle electromechanical delay
J = 1;                      % number of joints
M = 4;                      % number of muscles
S = 5;                      % number of muscle states
C = 4;                      % number of constraints of the model dynamics


%% load the experimental data
% initialize the matrix of optimization inputs
mus_act = zeros(sum(N), M);
torque = zeros(sum(N), J);
lmt = zeros(sum(N), M);
ma = zeros(sum(N), M);
hs = zeros(1, T);

for t = 1:T  % load trial data and get muscle parameters from the averaged gait data

    trial = trialNames(t);

    load(sprintf('D:/HuaweiWang/ReflexIdParaData/Subj%02d/Subj%02d_%s.mat', ...
        subj, subj, trial));

    % get muscle activation
    mus_act(sum(N(1:t-1))+1:sum(N(1:t)), :) = idParaData.mus_act(1:N(t), 1:M);
    
    % get joint torques
    torque(sum(N(1:t-1))+1:sum(N(1:t)), :) = idParaData.torque(1:N(t), 1:J);

    % get muscle length and moment arms
    lmt(sum(N(1:t-1))+1:sum(N(1:t)), :) = idParaData.lmt(1:N(t), 1:M);
    ma(sum(N(1:t-1))+1:sum(N(1:t)), :) = idParaData.ma(1:N(t), 1:M);
    
    % get the time interval
    hs(t) = idParaData.hs;
    
end

par0 = idParaData.mus_par0; 
P = length(par0);


%% generate initial guesses
% initialize the optimizing parameters
mus_a = mus_act(1:sum(N), :);
mus_da = zeros(size(mus_a));
for t = 1:T
    if t== 1
        mus_da(2:N(t), :) = (mus_act(2:N(t), :) - mus_act(1:N(t)-1, :))./hs(t);
    else
        mus_da(sum(N(1:t-1))+2:sum(N(1:t)), :) = ...
            (mus_act(sum(N(1:t-1))+2:sum(N(1:t)), :) ...
            - mus_act(sum(N(1:t-1))+1:sum(N(1:t))-1, :))./hs(t);
    end
end

lce = 0.05 + 0.02*rand(size(lmt(1:sum(N), :)));
dlce = -0.05 + 0.1*rand(size(lmt(1:sum(N), :)));

x0_1 = [mus_a(1:sum(N), :), mus_da(1:sum(N), :),...
        lce(1:sum(N), :), dlce(1:sum(N), :), mus_act(1:sum(N), :)];

x0_2 = reshape(x0_1', [1, sum(N)*M*S]);

par = zeros(1, P);  % randomize muscle parameters
par(1:M) = (0.9 + 0.2*rand(1, M)).*par0(1:M);  % lce_opt
par(M+1:2*M) = (0.9 + 0.2*rand(1, M)).*par0(M + 1:2*M);  % lt_slack
par(2*M+1:3*M) = par0(M + 1:2*M);  % theta0
par(3*M+1:4*M) = (0.9 + 0.2*rand(1, M)).*par0(3*M+1:4*M);  % Fmax

x = [x0_2, par];

Jac = jacobian_MPO(x, M, S, C, N, P, lmt, t_em, hs);
[row, col] = jacobianStructure_MPO(M, S, C, N, P, lmt, t_em, hs);

rCons = M*C*sum(N - ceil(t_em./hs) - 1) ...
            + M*(C-1)*sum(ceil(t_em./hs));

Jac_equ = zeros(rCons, sum(N)*M*S + P);

for i = 1:length(Jac)
    Jac_equ(row(i), col(i)) = Jac(i);
end

% finite differentiation
delta = 1e-6;


% differentiation of muscle fiber lengths lce
for ia = 1:length(x)
    
   % get values with upper change
   delta_x = x(ia)*delta;
   
   x(ia) = x(ia) + delta_x;
   cons_up = constraints_MPO(x, M, S, C, N, P, lmt, t_em, hs);
   
   % get values with upper change
   x(ia) = x(ia) - 2*delta_x;
   cons_dw = constraints_MPO(x, M, S, C, N, P, lmt, t_em, hs);

   % change back to the original value
   x(ia) = x(ia) + delta_x;

   % calculate the finite differentiation
   Jac_fd(:, ia) = (cons_up - cons_dw)/(2*delta_x);

end

% check the drivative differences
tolerance = 1e-5;

% difference check of df_da
errorid_Jac = diffEvaluate(Jac_equ, Jac_fd, tolerance);

if ~isempty(errorid_Jac)
   fprintf('Differentiations in Jac beyond thresholds\n')
end