function [errorid, realDiff, relativeDiff] = diffEvaluate(df_dx, df_dx_fd, tolerance)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate differences between two derivative matrix
    %
    % By: Huawei Wang
    % Date: 12/06/2021
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    realDiff = abs(df_dx - df_dx_fd);  % real difference errors
    relativeDiff = zeros(size(realDiff));
    
    nonzeroid = find(df_dx);  % nonzero element id in df_dx
    % relative difference errors
    relativeDiff(nonzeroid) = abs(realDiff(nonzeroid)./df_dx(nonzeroid));  
    
    % find the indexes of errors that larger than the tolerance
    errorid_real = find(realDiff > tolerance);
    errorid_diff = find(relativeDiff > tolerance);
    
    errorid = intersect(errorid_real, errorid_diff);

end