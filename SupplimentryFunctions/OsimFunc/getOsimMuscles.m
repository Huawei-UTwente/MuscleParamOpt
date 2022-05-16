function osimMuscleList = getOsimMuscles(OsimModel, MuscleNames)

    osimMuscleList = [];
    osimMuscleOrder = [];

    % get frame list
    muscleList = OsimModel.getMuscleList();

    % select frame based on the body name
    iter = muscleList.begin();

    while iter ~= muscleList.end() % run through all frames until get the same body name
        if ismember(string(iter.getName()), MuscleNames)
            osimMuscleList = [osimMuscleList, iter.deref()];
            osimMuscleOrder = [osimMuscleOrder, find(contains(MuscleNames, string(iter.getName())))];
        end
        
        if length(osimMuscleOrder) == length(MuscleNames)
            break;
        else
            iter.next();
        end
    end
    
    osimMuscleList = osimMuscleList(osimMuscleOrder);
    
end