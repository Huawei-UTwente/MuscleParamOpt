%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative test of the contaction dynamics
%
% By: Huawei Wang
% Date: 12/06/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

% initlize parameters
T = 1;  % number of data trials
N = [100]; % number of data node at each data trial
J = 1;  % number of joints
M = 4;  % number of muscles
S = 5;  % number of muscle states
C = 5;  % number of constraints of the model dynamics
tPhase = max(N); % number of totoal phase of the reflex controller
phase = {1:N(1)}';  % phase of the reflex controller
nPar_rf = 2*tPhase*M*M + 3*M;

t_em = 0.0 + zeros(1, T);
t_rf = 0.0 + zeros(1, T);

trialNames = ["walk_09", "walk_18", "walk_27", "walk_36", "walk_45", "walk_54", "run_63", "run_81", "run_99"];
subjMass = [69.3, 97.8, 58.4, 75.2, 65.1, 56.0, 80.8, 61.5, 81, 61.6, 62.4, 69];
muscleNames = ["soleus_r", "lat_gas_r", "med_gas_r", "tib_ant_r"];
coordinateNames = ["ankle_angle_r"];


for subj = 6
    
    % get subject mass
    mass = subjMass(subj);
    
    % load opensim model
    OsimModelFile = sprintf('../../PortableSystem_DataAnalysis/Processed_data/Subj%02d/OS/UTmodel/gait2392_simbody_subj%02d_scaled_1.osim', subj, subj);
    
    save_path = sprintf('ReflexControllerOptResults/Subj%02d', subj);
    
    % get muscle parameters
    mus_save_path = sprintf('MuscleParameterOptResults/Subj%02d', subj);
    par_mus = load(sprintf('%s/mus_par.mat', mus_save_path));
    par_mus = par_mus.mus_par;
    
    par_mus(3*M + 1:4*M) = par_mus(3*M + 1:4*M)/mass;
    
    % initialize the matrix of optimization inputs
    mus_act = zeros(sum(N), M);
    moment = zeros(sum(N), J);
    lmt = zeros(sum(N), M);
    ma = zeros(sum(N), M);
    hs = zeros(1, T);
    
    for t = 1:T  % load trial data and get muscle parameters from the averaged gait data
        
        trial = trialNames(t);
        
        processedData = importdata(sprintf('../../PortableSystem_DataAnalysis/Processed_data/Subj%02d/Subj%02d_%s.mat', subj, subj, trial));

        % get muscle activation
        muscleIndex = [];
        for musclename = muscleNames(1:M)
            muscleIndex = [muscleIndex, find(processedData.EMG.DataLabel == musclename)];
        end
        if t == 1
            mus_act(1:N(t), :) = processedData.Resample.Sych.Average.EMG.ave_r(1:N(t), muscleIndex+1);
        else
            mus_act(sum(N(1:t-1))+1:sum(N(1:t)), :) = processedData.Resample.Sych.Average.EMG.ave_r(1:N(t), muscleIndex+1);
        end
        % get joint torques
        jointIndex = [];
        for coordinatename = coordinateNames(1:J)
            for coori = 1:length(processedData.Resample.Sych.IDTrqDataLabel)
                if contains(processedData.Resample.Sych.IDTrqDataLabel{coori}, coordinatename)
                    jointIndex = [jointIndex, coori];
                end
            end
        end
        if t == 1
            moment(1:N(t), :) = processedData.Resample.Sych.Average.IDTrqData.ave_r(1:N(t), jointIndex)/mass;
        else
            moment(sum(N(1:t-1))+1:sum(N(1:t)), :) = processedData.Resample.Sych.Average.IDTrqData.ave_r(1:N(t), jointIndex)/mass;
        end
        
        % get muscle length and moment arms
        ikData.data = processedData.Resample.Sych.Average.IKAngData.ave_r(1:N(t), :);
        ikData.colheaders = processedData.Resample.Sych.IKAngDataLabel;
        if t == 1
            [lmt(1:N(t), :), ma(1:N(t), :)] = ...
            getOsimMuscleLengthMA(OsimModelFile, ikData, muscleNames(1:M), coordinateNames(1:J));
        else
            [lmt(sum(N(1:t-1))+1:sum(N(1:t)), :), ma(sum(N(1:t-1))+1:sum(N(1:t)), :)] = ...
            getOsimMuscleLengthMA(OsimModelFile, ikData, muscleNames(1:M), coordinateNames(1:J));
        end
        % get the time interval
        hs(t) = mean(processedData.Resample.Sych.Average.hsMatrix_right(:, 2) - processedData.Resample.Sych.Average.hsMatrix_right(:, 1))/100/100;
        
    end

    %% generate initial guesses
    
    opt = 1;
    
    load(sprintf('%s/optimization_res%02d.mat', save_path, opt));
    
end

x = [reshape(states', 1, N*M*S), parameters];

cons = constraints_RPO3(x, M, S, C, N, nPar_rf, lmt, par_mus, t_em, t_rf, hs, phase, tPhase);