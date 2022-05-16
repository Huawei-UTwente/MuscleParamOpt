
function [OsimModel, OsimState] = updOsimStateVariableValue_id(OsimModel,...
                                OsimState, kinematicsData, row)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Update state variables inside the OsimModel based on given kinematics
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get the Opensim model state names
    OsimStateNames = OsimModel.getStateVariableNames();

    for i = 0:OsimStateNames.getSize()-1  % run over all States
       stateName = string(OsimStateNames.get(i));  % get current state name
       bshIndex = strfind(stateName, '/');  % get the state type

       % if the state type is 'value' or 'velocity', find the location of
       % corresponding variables in the kinematics data
       if extractBetween(stateName, bshIndex(end)+1, strlength(stateName)) == "value"

           % find the location of the state name inside the variable list
           varIndexKin = find(kinematicsData.varNames == extractBetween(stateName, bshIndex(end-1)+1, bshIndex(end)-1));

           if  ~isempty(varIndexKin) % if the location is not empty
               
               % set the given value of this state
               OsimModel.setStateVariableValue(OsimState, OsimStateNames.get(i), kinematicsData.varValues(row, varIndexKin));
               
           end
       end

    end

end
