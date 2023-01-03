%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-processing of the gait data
%
% By: Huawei Wang
% Date: 04/05/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

dataFile = '../D2.4/SecondRecording/Processed/Processed_lift_heel_rise.mat';
frameSelect = 16001:3:16600;
osimModelFile = '../D2.4/SecondRecording/Processed/OS/UTmodel/gait2392_Demo_scaled.osim';
muscleNames = ["soleus_r", "lat_gas_r", "med_gas_r", "tib_ant_r"];
coordNames = ["ankle_angle_r"];
mass = 70;
savingPath = 'D2.4';

idDataPreparison(dataFile, frameSelect, osimModelFile, ....
    muscleNames, coordNames, mass, savingPath)

function idDataPreparison(dataFile, frameSelect, osimModelFile, ....
    muscleNames, coordNames, mass, savingPath)

    if ~exist(savingPath, 'dir')
        mkdir(savingPath)
    end

    % general parameters
    M = length(muscleNames);
    J = length(coordNames);

    % load muscle parameters
    

    % load propcessed experimental data
    processedData = importdata(dataFile);

    % get muscle activation
    muscleIndex = [];
    for musclename = muscleNames
        muscleIndex = [muscleIndex, find(processedData.Resample.Sych.NormEMGLabel == musclename)];
    end
        mus_act = processedData.Resample.Sych.NormEMGData(frameSelect, muscleIndex);
    
    % get joint torques
    jointIndex = [];
    for coordname = coordNames
        for coori = 1:length(processedData.Resample.Sych.IDTrqDataLabel_Portable)
            if contains(processedData.Resample.Sych.IDTrqDataLabel_Portable{coori}, coordname)
                jointIndex = [jointIndex, coori];
            end
        end
    end

    torque = processedData.Resample.Sych.IDTrqData_Portable(frameSelect, jointIndex)/mass;

    % get muscle length and moment arms
    ikData.data = processedData.Resample.Sych.IMUAngData(frameSelect, :);
    ikData.colheaders = processedData.Resample.Sych.IMUAngDataLabel;
    
    [lmt, ma] = getOsimMuscleLengthMA(osimModelFile, ikData, ...
        muscleNames(1:M), coordNames);
    
    % get the time interval
    hs = processedData.Resample.Sych.IMUAngData(frameSelect(2), 1) - processedData.Resample.Sych.IMUAngData(frameSelect(1), 1);

    % get ankle joint angle index
    for coori = 1:length(processedData.Resample.Sych.IMUAngDataLabel)
        if contains(processedData.Resample.Sych.IMUAngDataLabel{coori}, 'ankle_angle_r')
            ankle_ik_id = coori;
        end
    end

    [fiberOpt, tendonSlack, pennaAng, forceMax] =...
    getOsimMuscleParameter(osimModelFile, muscleNames);
    
    idParaData.mus_act = mus_act;
    idParaData.mus_par = [fiberOpt, tendonSlack, pennaAng, forceMax];
    idParaData.torque = torque;
    idParaData.lmt = lmt;
    idParaData.ma = ma;
    idParaData.hs = hs;
    idParaData.mass = mass;
    idParaData.mus_names = muscleNames;
    idParaData.coord_names = coordNames;
    
    dataNames = strsplit(dataFile, '/');
    save(sprintf('%s/%s', savingPath, dataNames{end}), 'idParaData');

end
