function [fiberOpt, tendonSlack, pennaAng, forceMax] =...
    getOsimMuscleParameter(OsimModelFile, muscleNames)
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

    % select left calcn frame
    Muscles = getOsimMuscles(OsimModel, muscleNames);

    % lengths of output data rows and columns
    col = length(muscleNames);

    % initilize parameter result vectors
    fiberOpt = zeros(1, col);
    tendonSlack = zeros(1, col);
    forceMax = zeros(1, col);
    pennaAng = zeros(1, col);
    
    % extract muscle parameters
    for c = 1:col
        fiberOpt(c) = Muscles(c).get_optimal_fiber_length();
        tendonSlack(c) = Muscles(c).get_tendon_slack_length();
        forceMax(c) = Muscles(c).get_max_isometric_force();
        pennaAng(c) = Muscles(c).get_pennation_angle_at_optimal();
    end
        
end

