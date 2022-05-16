function osimJointList = getOsimJoints(OsimModel, jointNames)

    osimJointList = [];
    osimJointOrder = [];

    % get frame list
    jointList = OsimModel.getJointList();

    % select frame based on the body name
    iter = jointList.begin();

    while iter ~= jointList.end()   % run through all frames until get the same body name
        if ismember(string(iter.getName()), jointNames)
            osimJointList = [osimJointList, iter.deref()];
            osimJointOrder = [osimJointOrder, find(contains(jointNames, string(iter.getName())))];
        end
        
        if length(osimJointOrder) == length(jointNames)
            break;
        else
            iter.next();
        end
    end
    
    osimJointList = osimJointList(osimJointOrder);
    
end

