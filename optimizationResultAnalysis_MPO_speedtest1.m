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

homeDataPath = 'D2.4';
trialNamesOpt = ["Processed_walk"];
muscleNames = ["soleus_l","lat_gas_l","med_gas_l","tib_ant_l"];
coordinateNames = ["ankle_angle_l"];
homeSavingPath = 'D2.4/walk';
mass = 69;

optimizationResultAnalysis_MPO_func(homeDataPath, trialNamesOpt, ...
                                        mass, muscleNames, coordinateNames,...
                                        homeSavingPath)


function optimizationResultAnalysis_MPO_func(homeDataPath, trialNamesOpt, ...
                                        mass, muscleNames, coordinateNames,...
                                        homeSavingPath)

    T = length(trialNamesOpt);  % number of data trails involved in the optimization, each data trial
            % could have different data nodes.
    
    t_em = 0.1 + zeros(1, T);  % electrialmechaincal delay inside muscle physiologies, around
                               % 100ms. Could be different among data trials
    J = length(coordinateNames);  % number of joint in each side of leg
    M = length(muscleNames);  % number of muscle in each side of leg
    S = 5;  % number of states of each muscle model
    C = 4;  % number of constraints of each muscle model

    % initialize the matrix of optimization inputs
%     mus_act = zeros(sum(N), M);
%     torque = zeros(sum(N), J);
%     lmt = zeros(sum(N), M);
%     d = zeros(sum(N), M);
    hs = zeros(1, T);

    for t = 1:T  % load trial data and get muscle parameters from the averaged gait data

        trial = trialNamesOpt(t);

        musPar = load(sprintf('%s/%s.mat', ...
            homeDataPath, trial));

        % get muscle activation
        mus_act = musPar.idParaData.mus_act(:, 1:M);

        % get joint torques
        torque = musPar.idParaData.torque(:, 1:J);

        % get muscle length and moment arms
%         lmt = musPar.idParaData.lmt(:, 1:M);
%         d = musPar.idParaData.ma(:, 1:M);

        % get the time interval
        hs(t) = musPar.idParaData.hs;

    end
    
    N = length(torque);

    % initial muscle parameters from the scaled Osim model
    mus_par0 = musPar.idParaData.mus_par; 
    P = length(mus_par0);
    
    % optimization saving path
    mkdir(homeSavingPath);

    % initial muscle parameters from the scaled Osim model
    % mus_par0 = musPar.idParaData.mus_par0;
    lce_opt0 = mus_par0(1:M);
    lt_slack0 = mus_par0(M + 1:2*M);
    theta0 = mus_par0(2*M + 1:3*M);
    Fmax0 = mus_par0(3*M + 1:4*M);
    
    save_path = homeSavingPath;

    % weights of the terms of the objective function
    W1 = 50;  % weight of joint torque fits
    W2 = 100;  % weight of muscle activation fits
    W3 = 10;  % weight of muscle activation smoothness
    W4 = 10;  % weight of muscle force smoothness
    W5 = 10; % weight of diversity of the optimizing parameters

    Range = 0.25;

    %% load optimized results

    for opt = 1:86

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
    % ylim([0, 5])

    st = 1;
    ed = length(Obj_sort);
    subplot(3, 1, 3)
    plot([1, 14], [1, 1],...
        '-', 'linewidth', 1, 'Color',[0 0 0]+0.75)
    hold on
    boxplot([Lce_opt_sort(st:ed, :)./lce_opt0, Lt_slack_sort(st:ed, :)./lt_slack0, ...
            Theta0_sort(st:ed, :)./theta0, Fmax_sort(st:ed, :)./Fmax0],...
        'Labels',{'Lce_opt1', 'Lce_opt2', 'Lce_opt3', 'Lce_opt4', ...
        'Lt_slack1', 'Lt_slack2', 'Lt_slack3', 'Lt_slack4', ...
        'Theta_opt1', 'Theta_opt2', 'Theta_opt3', 'Theta_opt4',...
        'Fmax1', 'Fmax2', 'Fmax3', 'Fmax4'})
    ylim([1-2*Range, 1+4*Range])
    hold off
    title('Optimized Parameters')
    ylabel('Normalized range %')
    xtickangle(angle)

    st = 1;
    ed = 1; %min(ceil(length(Obj_sort)/2), 3); %length(Obj_sort);
    subplot(3, 1, 2)
    plot([1, 14], [1, 1],...
        '-', 'linewidth', 1, 'Color',[0 0 0]+0.75)
    hold on
    plot([Lce_opt_sort(st:ed, :)./lce_opt0, Lt_slack_sort(st:ed, :)./lt_slack0,...
          Theta0_sort(st:ed, :)./theta0, Fmax_sort(st:ed, :)./Fmax0]', '*')
    ylim([1-2*Range, 1+4*Range])

    hold off
    title('Optimized Parameters in top 3 optimizations')
    ylabel('Normalized range %')
    xtickangle(angle)
    savefig(fig4, sprintf('%s/objectiveFunction_OptimizedParameters.fig', save_path))

    % save the calibrated (best) muscle paraemters
    mus_par = [Lce_opt_sort(1, :), Lt_slack_sort(1, :), Theta0_sort(1, :), Fmax_sort(1, :)];
    save(sprintf('%s/mus_par.mat', save_path), 'mus_par');

    % plot joint torque comparison
    plot_ind = succ_index(sort_index);

    color_vec = ['r', 'b', 'g', 'm', 'c', 'y'];

    color_ind = 1;

    h1 = figure(1);
    joint_names = ["ANKLE MOMENT Nm"];
    
    best_res = load(sprintf('%s/optimization_res%02d.mat', save_path, plot_ind(1)));
    
    for opt = plot_ind(st:ed)

        saving_names = sprintf('%s/optimization_res%02d.mat', save_path, opt);
        res = load(saving_names);

        subj_h1 = ceil(sqrt(J));

        for j = 1:J
            subplot(subj_h1, subj_h1, j)
            for t = 1:T
                plot(1:N(t), torque(sum(N(1:t-1))+1:sum(N(1:t)), :)*mass,...
                    'k-', 'linewidth', 2.5)
                hold on
                plot(1:N(t), res.mom_res(sum(N(1:t-1))+1:sum(N(1:t)), :)*mass,...
                    '--', 'linewidth', 2.5, 'color', color_vec(color_ind))
                hold on

                if t == T
                    title(joint_names(j)) 
                    legend("experimental Joint Torques", ...
                           "optimized Joint Torques")
                end
                
                r_Tor = corrcoef(torque(sum(N(1:t-1))+1:sum(N(1:t)), :), ...
                    res.mom_res(sum(N(1:t-1))+1:sum(N(1:t)), :));
                
                best_res.coeff.torque(t) = r_Tor(1, 2);
                best_res.rms.torque(t) = rms(torque(sum(N(1:t-1))+1:sum(N(1:t)), :) ...
                    - res.mom_res(sum(N(1:t-1))+1:sum(N(1:t)), :));
                
            end
        end
        color_ind = color_ind + 1;
    end
    savefig(h1, sprintf('%s/JointTorqueFits.fig', save_path))
    
    close(h1)

    h2 = figure(2);
    color_ind = 1;

    for opt = plot_ind(st:ed)

        saving_names = sprintf('%s/optimization_res%02d.mat', save_path, opt);
        res = load(saving_names);

        for t = 1:T
            for m = 1:M
                subplot(T, M, (t-1)*M + m)
                plot(1:N(t), res.force_res(sum(N(1:t-1))+1:sum(N(1:t)), m),...
                    '--', 'linewidth', 2.5, 'color', color_vec(color_ind))
                hold on

                if t == 1
                    title(muscleNames(m))
                end

                if m == 1
                    trail_name = trialNamesOpt(t);
                    ylabel(trail_name)
                end
            end
            
        end
        color_ind = color_ind + 1;
        sgtitle('Muscle Forces')
    end
    savefig(h2, sprintf('%s/MuscleForces.fig', save_path))
    close(h2)
    
    %%    
    h3 = figure(3);
    color_ind = 1;

    for opt = plot_ind(st:ed)

        saving_names = sprintf('%s/optimization_res%02d.mat', save_path, opt);
        res = load(saving_names);

        for t = 1:T
            for m = 1:M
                subplot(T, M, (t-1)*M + m)
                plot(1:N(t), mus_act(sum(N(1:t-1))+1:sum(N(1:t)), m),...
                    'k-', 'linewidth', 2.5)
                hold on
                plot(1:N(t), res.states(sum(N(1:t-1))+1:sum(N(1:t)), 4*M + m),... %/max_act(m)*0.3,...
                    '--', 'linewidth', 2.5, 'color', color_vec(color_ind))
                
                ylim([0, 0.7])

                if t == 1
                    title(muscleNames(m))
                end

                if m == 1
                    trail_name = trialNamesOpt(t);
                    ylabel(trail_name)
                end
                
                r_Mus = corrcoef(mus_act(sum(N(1:t-1))+1:sum(N(1:t)), m), ...
                    res.states(sum(N(1:t-1))+1:sum(N(1:t)), 4*M + m));
                
                best_res.coeff.activation(t, m) = r_Mus(1, 2);
                best_res.rms.activation(t, m) = rms(mus_act(sum(N(1:t-1))+1:sum(N(1:t)), m) ...
                    - res.states(sum(N(1:t-1))+1:sum(N(1:t)), 4*M + m));
                
            end
            sgtitle('Muscle Activation')
        end

        color_ind = color_ind + 1;
    end
    savefig(h3, sprintf('%s/MuscleActivation.fig', save_path))
    close(h3)
    
    h4 = figure();
    boxplot(Time_res_sort);
    savefig(h4, sprintf('%s/ComputingTime.fig', save_path))
    close(h4)
    
    save(sprintf('%s/best_res.mat',save_path), 'best_res');
    
end