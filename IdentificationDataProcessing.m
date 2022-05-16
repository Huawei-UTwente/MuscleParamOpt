%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-processing of the gait data
%
% By: Huawei Wang
% Date: 04/05/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

dataFile = 'ExperimentalData/walk_36.mat';
osimModelFile = 'ExperimentalData/gait2392.osim';
MuscleParaFile = 'ExperimentalData/par_mus.mat';
muscleNames = ["soleus_r", "lat_gas_r", "med_gas_r", "tib_ant_r"];
coordNames = ["ankle_angle_r"];
mass = 58.4;
savingPath = 'idParaData';

idDataPreparison(dataFile, osimModelFile, MuscleParaFile, ....
    muscleNames, coordNames, mass, savingPath)

function idDataPreparison(dataFile, osimModelFile, MuscleParaFile, ....
    muscleNames, coordNames, mass, savingPath)

    if ~exist(savingPath, 'dir')
        mkdir(savingPath)
    end

    % general parameters
    M = length(muscleNames);
    J = length(coordNames);

    % load muscle parameters
    load(MuscleParaFile);
    par_mus(3*M + 1:4*M) = par_mus(3*M + 1:4*M)/mass;  % normalize Fmax with mass

    % load propcessed experimental data
    processedData = importdata(dataFile);

    % get muscle activation
    muscleIndex = [];
    for musclename = muscleNames
        muscleIndex = [muscleIndex, find(processedData.EMG.DataLabel == musclename)];
    end
        mus_act = processedData.Resample.Sych.Average.EMG.ave_r(:, muscleIndex+1);
    
    % get joint torques
    jointIndex = [];
    for coordname = coordNames
        for coori = 1:length(processedData.Resample.Sych.IDTrqDataLabel)
            if contains(processedData.Resample.Sych.IDTrqDataLabel{coori}, coordname)
                jointIndex = [jointIndex, coori];
            end
        end
    end

    torque = processedData.Resample.Sych.Average.IDTrqData.ave_r(:, jointIndex)/mass;

    % get muscle length and moment arms
    ikData.data = processedData.Resample.Sych.Average.IKAngData.ave_r;
    ikData.colheaders = processedData.Resample.Sych.IKAngDataLabel;
    
    [lmt, ma] = getOsimMuscleLengthMA(osimModelFile, ikData, ...
        muscleNames(1:M), coordNames);
    
    % get the time interval
    hs = mean(processedData.Resample.Sych.Average.hsMatrix_right(:, 2) - processedData.Resample.Sych.Average.hsMatrix_right(:, 1))/100/100;

    % get ankle joint angle index
    for coori = 1:length(processedData.Resample.Sych.IKAngDataLabel)
        if contains(processedData.Resample.Sych.IKAngDataLabel{coori}, 'ankle_angle_r')
            ankle_ik_id = coori;
        end
    end
    
    % get Fy index
    for coori = 1:length(processedData.Resample.Sych.ForcePlateGRFDataLabel)
        if contains(processedData.Resample.Sych.ForcePlateGRFDataLabel{coori}, 'ground_force_vy')
            fy_id = coori;
        end
    end

    phase = getPhase(processedData.Resample.Sych.Average.IKAngData.ave_r(:, ankle_ik_id), ...
        processedData.Resample.Sych.Average.ForcePlateGRFData.ave_r(:, fy_id)/mass);
    
    idParaData.mus_act = mus_act;
    idParaData.mus_par = par_mus;
    idParaData.torque = torque;
    idParaData.lmt = lmt;
    idParaData.ma = ma;
    idParaData.hs = hs;
    idParaData.mass = mass;
    idParaData.phase = phase;
    idParaData.mus_names = muscleNames;
    idParaData.coord_names = coordNames;
    
    dataNames = strsplit(dataFile, '/');
    save(sprintf('%s/%s', savingPath, dataNames{end}), 'idParaData');

end
