function osimCoordinateList = getOsimCoordinates(OsimModel, coordinateNames)

    osimCoordinateList = [];
    osimCoordinateOrder = [];

    % get coordinate list
    coordinateList = OsimModel.getCoordinateSet();
    coordinateNum = OsimModel.getNumCoordinates();
    
    for iter = 0:coordinateNum   % run through all coordinate until get the same name
        if ismember(string(coordinateList.get(iter)), coordinateNames)
            osimCoordinateList = [osimCoordinateList, coordinateList.get(iter)];
            osimCoordinateOrder = [osimCoordinateOrder,...
                find(contains(coordinateNames, string(coordinateList.get(iter))))];
        end
        
        if length(osimCoordinateOrder) == length(coordinateNames)
            break;
        end
    end
    osimCoordinateList = osimCoordinateList(osimCoordinateOrder);
    
end

