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

save_path = 'MuscleParameterOptResults/Subj06';

for subj = 06
    
    musPar = [];
    
    musPar.info.subjLevel.subjNum = "subject number in the UT walking and running experiment, int";
    musPar.info.subjLevel.mass = "subject mass, kg";
    musPar.info.subjLevel.musName = "list of muscle names examed, string";
    musPar.info.subjLevel.jntName = "list of joint names examed, string";
    
    musPar.info.trialLevel.lmt = "muscle tendon unit lengths of examed muscles corresponding to the trial kinematics, m";
    musPar.info.trialLevel.ma = "moment arms of examed muscles corresponding to the trial kinematics, in m";
    musPar.info.trialLevel.h = "time interval of this data trial (different between trials), second";

    musPar.subjNum = subj;
    musPar.mass = subjMass(subj);
    musPar.musName = muscleNames;
    musPar.jntName = coordinateNames;
    
    % load opensim model
    OsimModelFile = sprintf('../../PortableSystem_DataAnalysis/Processed_data/Subj%02d/OS/UTmodel/gait2392_simbody_subj%02d_scaled_1.osim', subj, subj);
    
%     % get muscle parameters
%     [lce_opt0, lt_slack0, theta0, Fmax0] = ...
%             getOsimMuscleParameter(OsimModelFile, muscleNames(1:M));
    
    % initialize the matrix of optimization inputs
    mus_act = zeros(N, M);
    moment = zeros(N, J);
    lmt = zeros(N, M);
    ma = zeros(N, M);
    
    for t = 1:T  % load trial data and get muscle parameters from the averaged gait data
        
        trial = trialNames(t);
        
        processedData = importdata(sprintf('../../PortableSystem_DataAnalysis/Processed_data/Subj%02d/Subj%02d_%s.mat', subj, subj, trial));

        % get muscle length and moment arms
        ikData.data = processedData.Resample.Sych.Average.IKAngData.ave_r(delay_rm+1:N+delay_rm, :);
        ikData.colheaders = processedData.Resample.Sych.IKAngDataLabel;
        [lmt, ma] = ...
            getOsimMuscleLengthMA(OsimModelFile, ikData, muscleNames(1:M), coordinateNames(1:J));
        
        % get the time interval
        h = mean(processedData.Resample.Sych.Average.hsMatrix_right(:, 2) - processedData.Resample.Sych.Average.hsMatrix_right(:, 1))/100/100;
        
        musPar.(trial).lmt = lmt;
        musPar.(trial).ma = ma;
        musPar.(trial).h = h;

    end
end

%% load optimized results

for opt = 1:100
    
    saving_names = sprintf('%s/optimization_res%02d.mat', save_path, opt);
    res = load(saving_names);
    Obj_res(opt) = res.obj;
    Status_res(opt) = res.status;
    
end

%% plot optimized parameters
plot_sign = 0;
succ_index = find(Status_res == plot_sign);
[Obj_sort, sort_index] = sort(Obj_res(succ_index));

best_idx_sort = succ_index(sort_index);


saving_names = sprintf('%s/optimization_res%02d.mat', save_path, best_idx_sort(1));
res = load(saving_names);

musPar.optFiber = res.parameters(1:4);
musPar.optSlack = res.parameters(5:8);
musPar.optPenna = res.parameters(9:12);
musPar.optFmax = res.parameters(13:16)*musPar.mass;

musPar.info.trialLevel.optFiber = "Optimal fiber lengths, m";
musPar.info.trialLevel.optSlack = "tendon slack lengths, m";
musPar.info.trialLevel.optPenna = "pennation angle at the optimal fiber length, rad";
musPar.info.trialLevel.optFmax = "Fmax of examed muscles, N";

for t = 1:T
    trial = trialNames(t);
    
    musPar.(trial).sti = res.activations(:, (t-1)*M+1:t*M);
    musPar.(trial).act = res.states(:, (t-1)*M+1:t*M);
    musPar.(trial).dact = res.states(:, T*M + (t-1)*M+1:T*M+ t*M);
    musPar.(trial).lce = res.states(:, 2*T*M + (t-1)*M+1:2*T*M + t*M);
    musPar.(trial).dlce = res.states(:, 3*T*M + (t-1)*M+1:3*T*M+ t*M);
    
    musPar.(trial).musForce = res.force_res(:, (t-1)*M + 1:t*M);
    
end

musPar.info.trialLevel.sti = "stimulation of examed muscles, none";
musPar.info.trialLevel.act = "activations of examed muscles, none";
musPar.info.trialLevel.dact = "activation derivatives of examed muscles, none";
musPar.info.trialLevel.lce = "muscle fiber lengths of examed muscles, m";
musPar.info.trialLevel.dlce = "muscle fiber length derivatives of examed muscles, m/s";
musPar.info.trialLevel.musForce = "muscle forces of examed muscles, N";
musPar.info.trialLevel.trialName = "trial name indicate the locomotion type, the last two number indicate the speed, km/h";

save(sprintf('%s/musPar.mat', save_path), 'musPar')



