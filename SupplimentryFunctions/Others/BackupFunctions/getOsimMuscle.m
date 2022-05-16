function osimMuscle = getOsimMuscle(OsimModel, MuscleName)

    % get frame list
    muscleList = OsimModel.getMuscleList();

    % select frame based on the body name
    iter = muscleList.begin();

    while 1   % run through all frames until get the same body name
        if strcmp(iter.getName(),  MuscleName)
            osimMuscle = iter;
            break;
        end
        iter.next();
    end
   
end