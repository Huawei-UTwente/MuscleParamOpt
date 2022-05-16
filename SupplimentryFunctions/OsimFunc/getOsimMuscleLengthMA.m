function [lengthMTU, maAnkle] =...
    getOsimMuscleLengthMA(OsimModelFile, ikData, muscleNames, coordinateNames)
%  getOsimMuscleParameter? get muscle parameters from a OpenSim model,
%  especially the MTU lengths and moment arms will be calculated based on
%  kinematics data
%
% By: Huawei Wang
% Date: 12-04-2021

% load the Opensim libraries
import org.opensim.modeling.*

% load the OpenSim model
OsimModel = Model(OsimModelFile);

% % load the inverse kinematics results
% ikData = importdata(ikDataFile);

% get the data titles of kinematics data
varNames = string(ikData.colheaders);

% extract the values of coordinate variables (joint angles) and calculate
% corresponding velocities, using middle point differentiation 
varValues = ikData.data(:, 2:end);

% convert the kinematics joint angles from degree to radian. 4:6
% columns are the pelvis translation which does not need to convert.
% if contains(ikData.textdata{5, 1}, 'yes')  % if ik angles are in degrees
varValues(:, 1:3) = varValues(:, 1:3)*(pi/180);
varValues(:, 7:end) = varValues(:, 7:end)*(pi/180);

% combine kinematicsData together
kinematicsData.varNames = varNames(2:end);
kinematicsData.varValues = varValues;

% initilize the simulation model
OsimModelSim = initSystem(OsimModel);

% get the Opensim model states
OsimState = OsimModelSim.State;

% select left calcn frame
Muscles = getOsimMuscles(OsimModel, muscleNames);
Coordinates = getOsimCoordinates(OsimModel, coordinateNames);

% lengths of output data rows and columns
col = length(muscleNames);
row = length(varValues(:, 1));
hei = length(Coordinates);

lengthMTU = zeros(row, col);

if hei == 1
    maAnkle = zeros(row, col);
else
    maAnkle = zeros(row, col, hei);
end

for r = 1:row
    [OsimModel, OsimState] = updOsimStateVariableValue_id(OsimModel,...
                                OsimState, kinematicsData, r);
    for c = 1:col
            osimMuscle = Muscles(c);
            osimMuscle.computeEquilibrium(OsimState);
            lengthMTU(r, c) = osimMuscle.getLength(OsimState);
        if hei == 1
            maAnkle(r, c) = osimMuscle.computeMomentArm(OsimState, Coordinates);
        else
            for h = 1:hei
                maAnkle(r, c, h) = osimMuscle.computeMomentArm(OsimState, Coordinates(h));
            end
        end
        
    end
end

end

