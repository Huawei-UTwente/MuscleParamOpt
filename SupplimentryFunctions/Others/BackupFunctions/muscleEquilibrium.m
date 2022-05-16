function [musForce, musLength] = muscleEquilibrium(OsimModel, OsimState, ...
                                kinematicsData, activationData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Muscle Equilibrium Force by providing the kinematics and muscle
% activations. The current code does not consider the fiber velocities,
% therefore, only static movement is valid to use this code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % initilize the muscle force matrix
    musForce.varNames = activationData.varNames;
    musForce.varValues = zeros(size(activationData.varValues));
    
    musLength.varNames = activationData.varNames;
    musLength.varValues = zeros(size(activationData.varValues));
    
    % get the muscle model list from the given muscle activations
    osimMuscleList = [];
    for musName = activationData.varNames
        osimMuscleList = [osimMuscleList, getOsimMuscle(OsimModel, musName)];
    end
    
    % update the model states frame by frame, and get muscle forces
    if size(kinematicsData.varValues, 1) ~= size(activationData.varValues, 1)
        error('the lengths between kinematics and activations are not the same \n')
        
    else
        
        for row = 1:size(kinematicsData.varValues, 1)
            
            % Update OpenSim model states (joint angles, muscle activations)
            [OsimModel, OsimState] = updOsimStateVariableValue(OsimModel, OsimState,...
                                kinematicsData, activationData, row);
                            
            % calculate the muscle lengths and forces through the
            % 'computeEquilibrium' function in OpenSim
            for col = 1:length(osimMuscleList)
                osimMuscle = osimMuscleList(col);
                osimMuscle.computeEquilibrium(OsimState);
                musLength.varValues(row, col) = osimMuscle.getLength(OsimState);
%                 % check if the muscle length is too low
%                 if osimMuscle.getLength(OsimState) < 0.15
%                     continue
%                 end
                musForce.varValues(row, col) = osimMuscle.getTendonForce(OsimState);
            end
        
        end
    end
                      
end