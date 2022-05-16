%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is to optimize the ankle joint reflex control with the
% experimental data.
%
% The experimental data consists with joint angles, joint torques, and
% muscle activations (EMG).
%
% The ankle joint reflex will general muscle activations for the ankle
% muscles based on the reflex feedback. Song's (2015) feedback loop is took
% as reference. 
%
% The ankle joint muscles including: Soleus; Tibiais Anterior;
% M/L Gastrocneius. 
%
% Optimization is based on the averaged half gait cycle's data. gradient
% based optimizer is used here to find the best reflex control gains that
% can explain the experimental data.
%
% By: Huawei Wang
% Date: Augst 1, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all

% define global parameters
global lmt ma moment mus_act L J N M P S C T W1 W2 W3 h h_s coordinateNames muscleNames;

delay_rm = 6;  % electrialmechaincal delay inside muscle physiologies, remove
               % 5 nodes of the end of muscle activation, and remove 
               % 5 nodes of the beginning of joint torques

J = 1;  % number of joint in each side of leg
L = 1;  % 1: only one side of leg;  2: two sides of legs
M = 4;  % number of muscle in each side of leg
T = 9;  % number of data trails involved in the optimization, so far each trail should have the same number of nodes

N = 100 - delay_rm;  % number of nodes in each trail of data
S = 4;  % number of states of each muscle model
C = 4;  % number of constraints of each muscle model

trialNames = ["walk_09", "walk_18", "walk_27", "walk_36", "walk_45", "walk_54", "run_63", "run_81", "run_99"];
subjMass = [69.3, 97.8, 58.4, 75.2, 65.1, 56.0, 80.8, 61.5, 81, 61.6, 62.4, 69];
muscleNames = ["soleus_r", "lat_gas_r", "med_gas_r", "tib_ant_r"];
coordinateNames = ["ankle_angle_r"];

for subj = 06
    
    % get subject mass
    mass = subjMass(subj);
    
    % load opensim model
    OsimModelFile = sprintf('D:/HuaweiWang/PortableSystem_DataAnalysis/Processed_data/Subj%02d/OS/UTmodel/gait2392_simbody_subj%02d_scaled_1.osim', subj, subj);
    
    % get muscle parameters
    [lce_opt0, lt_slack0, theta0, Fmax0] = ...
            getOsimMuscleParameter(OsimModelFile, muscleNames(1:M));
    
    % initialize the matrix of optimization inputs
    mus_act = zeros(N, M*L*T);
    moment = zeros(N, J*L*T);
    lmt = zeros(N, M*L*T);
    ma = zeros(N, M*L*T);
    h = zeros(1, T*L*M);
    h_s = zeros(1, T*L*M*S);
    
    for t = 1:T  % load trial data and get muscle parameters from the averaged gait data
        
        trial = trialNames(t);
        
        processedData = importdata(sprintf('D:/HuaweiWang/PortableSystem_DataAnalysis/Processed_data/Subj%02d/Subj%02d_%s.mat', subj, subj, trial));

        % get muscle activation
        muscleIndex = [];
        for musclename = muscleNames(1:M)
            muscleIndex = [muscleIndex, find(processedData.EMG.DataLabel == musclename)];
        end
        mus_act(1:N, (t-1)*M*L+1:t*M*L) = processedData.Resample.Sych.Average.EMG.ave_r(1:N, muscleIndex+1);
        
        % get joint torques
        jointIndex = [];
        for coordinatename = coordinateNames(1:J)
            for coori = 1:length(processedData.Resample.Sych.IDTrqDataLabel)
                if contains(processedData.Resample.Sych.IDTrqDataLabel{coori}, coordinatename)
                    jointIndex = [jointIndex, coori];
                end
            end
        end
        moment(1:N, (t-1)*J*L+1:t*J*L) = processedData.Resample.Sych.Average.IDTrqData.ave_r(delay_rm+1:N+delay_rm, jointIndex)/mass;
        
        % get muscle length and moment arms
        ikData.data = processedData.Resample.Sych.Average.IKAngData.ave_r(delay_rm+1:N+delay_rm, :);
        ikData.colheaders = processedData.Resample.Sych.IKAngDataLabel;
        [lmt(1:N, (t-1)*M*L+1:t*M*L), ma(1:N, (t-1)*M*L+1:t*M*L)] = ...
            getOsimMuscleLengthMA(OsimModelFile, ikData, muscleNames(1:M), coordinateNames(1:J));
        
        % get the time interval
        h((t-1)*M*L+1:t*M*L) = mean(processedData.Resample.Sych.Average.hsMatrix_right(:, 2) - processedData.Resample.Sych.Average.hsMatrix_right(:, 1))/100;
        h_s((t-1)*M*L*S+1:t*M*L*S) = mean(processedData.Resample.Sych.Average.hsMatrix_right(:, 2) - processedData.Resample.Sych.Average.hsMatrix_right(:, 1))/100;

    end

    
    save_path = sprintf('D:/HuaweiWang/MuscleOptResults/MuscleParameterOptResults/Subj%02d', subj);

end

% weights of the terms of the objective function
W1 = 50;  % weight of joint torque fits
W2 = 100;  % weight of muscle activation fits
W3 = 10;  % weight of muscle activation smoothness
W4 = 10;  % weight of muscle force smoothness
W5 = 10; % weight of diversity of the optimizing parameters

Range = 0.25;

% generate optimizing muscle parameters
muscle_par = [lce_opt0, lt_slack0, theta0, Fmax0/mass];
P = length(muscle_par);

%% load optimized results

for opt = 1:100
    
    saving_names = sprintf('%s/optimization_res%02d.mat', save_path, opt);
    res = load(saving_names);
    
    Lce_opt_res(opt, :) = res.parameters(1:M);
    Lt_slack_res(opt, :) = res.parameters(M + 1:2*M);
    Theta0_res(opt, :) = res.parameters(2*M + 1:3*M);
    Fmax_res(opt, :) = res.parameters(3*M + 1:4*M);
    
    Obj_res(opt) = res.obj;
    Status_res(opt) = res.status;
    Time_res(opt) = res.time;
end

%% plot optimized parameters
plot_sign = 0;
succ_index = find(Status_res == plot_sign);
[Obj_sort, sort_index] = sort(Obj_res(succ_index));

Lce_opt_sort = Lce_opt_res(succ_index(sort_index), :);
Lt_slack_sort = Lt_slack_res(succ_index(sort_index), :);
Theta0_sort = Theta0_res(succ_index(sort_index), :);
Fmax_sort = Fmax_res(succ_index(sort_index), :);

Time_res_sort = Time_res(succ_index(sort_index))';

fig4 = figure(4);

angle = -45;

subplot(3, 1, 1)
plot(Obj_sort, 'ro')
title('Objective Values')
xlabel('Optimization trails')
ylim([0, 5])

st = 1;
ed = length(Obj_sort);
subplot(3, 1, 2)
plot([1, 14], [1, 1],...
    '-', 'linewidth', 1, 'Color',[0 0 0]+0.75)
hold on
boxplot([Lce_opt_sort(st:ed, :)./lce_opt0, Lt_slack_sort(st:ed, :)./lt_slack0, ...
        Theta0_sort(st:ed, :)./theta0, Fmax_sort(st:ed, :)./Fmax0*mass],...
    'Labels',{'Lce_opt1', 'Lce_opt2', 'Lce_opt3', 'Lce_opt4', ...
    'Lt_slack1', 'Lt_slack2', 'Lt_slack3', 'Lt_slack4', ...
    'Theta_opt1', 'Theta_opt2', 'Theta_opt3', 'Theta_opt4',...
    'Fmax1', 'Fmax2', 'Fmax3', 'Fmax4'})
ylim([1-2*Range, 1+2*Range])
hold off
title('Optimized Parameters')
ylabel('Normalized range %')
xtickangle(angle)

st = 1;
ed = 1; %min(ceil(length(Obj_sort)/2), 3); %length(Obj_sort);
subplot(3, 1, 3)
plot([1, 14], [1, 1],...
    '-', 'linewidth', 1, 'Color',[0 0 0]+0.75)
hold on
plot([Lce_opt_sort(st:ed, :)./lce_opt0, Lt_slack_sort(st:ed, :)./lt_slack0,...
      Theta0_sort(st:ed, :)./theta0, Fmax_sort(st:ed, :)./Fmax0*mass]', '*')
ylim([1-2*Range, 1+2*Range])

hold off
title('Optimized Parameters in top 3 optimizations')
ylabel('Normalized range %')
xtickangle(angle)
savefig(fig4, sprintf('%s/objectiveFunction_OptimizedParameters.fig', save_path))

% save the calibrated (best) muscle paraemters
mus_par = [Lce_opt_sort(1, :), Lt_slack_sort(1, :), Theta0_sort(1, :), Fmax_sort(1, :)*mass];
save(sprintf('%s/mus_par.mat', save_path), 'mus_par');

% plot joint torque comparison
plot_ind = succ_index(sort_index);

color_vec = ['r', 'b', 'g', 'm', 'c', 'y'];

color_ind = 1;
    
h1 = figure(1);
joint_names = ["ANKLE MOMENT Nm"];
    
for opt = plot_ind(st:ed)
 
    saving_names = sprintf('%s/optimization_res%02d.mat', save_path, opt);
    res = load(saving_names);
    
    subj_h1 = ceil(sqrt(J));
    
    for j = 1:J
        subplot(subj_h1, subj_h1, j)
        for t = 1:T
            for l = 1:L
                plot((l-1)*N + 1:l*N, moment(:, (t-1)*L*J + (l-1)*J + j),...
                    'k-', 'linewidth', 2.5)
                hold on
                plot((l-1)*N + 1:l*N, res.mom_res(:, (t-1)*L*J + (l-1)*J + j),...
                    '--', 'linewidth', 2.5, 'color', color_vec(color_ind))
                hold on
                
                if l == 1
                    plot([(l-1)*N + 1, (l-1)*N + N], [0, 0],...
                    'k-', 'linewidth', 1.5)
                    hold on
                elseif l == 2
                    plot([(l-1)*N + 1, (l-1)*N + N], [0, 0],...
                    '-', 'linewidth', 1.5, 'color', [1, 1, 1]*0.75)
                    hold on
                end

                if t == 1
                    title(joint_names(j))
                end

                if t == T
                    xlabel('RIGHT - LEFT')
                end
                if j == 1
                    trail_name = trialNames(t);
                    ylabel(trail_name)
                end            
                if t == T
                    legend("experimental Joint Torques", ...
                           "optimized Joint Torques")
                end
            end
        end
    end
    
    color_ind = color_ind + 1;
    
end
savefig(h1, sprintf('%s/JointTorqueFits.fig', save_path))


h2 = figure(2);
color_ind = 1;
    
for opt = plot_ind(st:ed)
    
    saving_names = sprintf('%s/optimization_res%02d.mat', save_path, opt);
    res = load(saving_names);

    for t = 1:T
        for m = 1:M
            for l = 1:L
                subplot(T, M, (t-1)*M + m)
                plot((l-1)*N + 1:l*N, res.force_res(1:N, (t-1)*L*M + (l-1)*M + m),...
                    '--', 'linewidth', 2.5, 'color', color_vec(color_ind))
                hold on
                if l == 1
                    plot([(l-1)*N + 1, (l-1)*N + N], [0, 0],...
                    'k-', 'linewidth', 1.5)
                    hold on
                elseif l == 2
                    plot([(l-1)*N + 1, (l-1)*N + N], [0, 0],...
                    '-', 'linewidth', 1.5, 'color', [1, 1, 1]*0.75)
                    hold on
                end

                if t == 1
                    title(muscleNames(m))
                end

                if m == 1
                    trail_name = trialNames(t);
                    ylabel(trail_name)
                end
            end
        end
        sgtitle('Muscle Forces')
    end
    
    color_ind = color_ind + 1;
    
end
savefig(h2, sprintf('%s/MuscleForces.fig', save_path))
%%    
h3 = figure(3);
color_ind = 1;
    
for opt = plot_ind(st:ed)
    
    saving_names = sprintf('%s/optimization_res%02d.mat', save_path, opt);
    res = load(saving_names);

    for t = 1:T
        max_act = max([res.activations(1:N, (t-1)*L*M + 1: (t-1)*L*M + M);...
                           res.activations(1:N, (t-1)*L*M + M + 1: (t-1)*L*M + L*M)]);
        for m = 1:M
            for l = 1:L
                subplot(T, M, (t-1)*M + m)
                plot((l-1)*N + 1:l*N, mus_act(1:N, (t-1)*L*M + (l-1)*M + m),...
                    'k-', 'linewidth', 2.5)
                hold on
                plot((l-1)*N + 1:l*N, res.activations(1:N, (t-1)*L*M + (l-1)*M + m),... %/max_act(m)*0.3,...
                    '--', 'linewidth', 2.5, 'color', color_vec(color_ind))
                if l == 1
                    plot([(l-1)*N + 1, (l-1)*N + N], [0, 0],...
                    'k-', 'linewidth', 1.5)
                    hold on
                elseif l == 2
                    plot([(l-1)*N + 1, (l-1)*N + N], [0, 0],...
                    '-', 'linewidth', 1.5, 'color', [1, 1, 1]*0.75)
                    hold on
                end
                
                ylim([0, 0.7])

                if t == 1
                    title(muscleNames(m))
                end

                if m == 1
                    trail_name = trialNames(t);
                    ylabel(trail_name)
                end
            end
        end
        sgtitle('Muscle Activation')
    end
    
    color_ind = color_ind + 1;
end
savefig(h3, sprintf('%s/MuscleActivation.fig', save_path))