%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative test of the gausssian function
%
% By: Huawei Wang
% Date: 12/06/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

% initlize parameters
% b1 = 0.815;
% b2 = 1.055;
% b3 = 0.162;
% b4 = 0.063;

b1 = 0.433;
b2 = 0.717;
b3 = -0.030;
b4 = 0.200;


lce_nor = [0.8, 1.2];  % deactivation time constant

% calculate activation dynamics derivatives
dfact_dlce_nor = gaussianFunctionAct_diff(lce_nor, b1, b2, b3, b4);

dfact_dlce_nor_equ = diag(dfact_dlce_nor);

% finite differentiation
delta = 1e-4;

% differentiation of muscle activation a
for ia = 1:length(lce_nor)
    
   delta_a = lce_nor(ia)*delta;
   
   % get values with upper change of a
   lce_nor(ia) = lce_nor(ia) + delta_a;
   fact_up = gaussianFunctionAct(lce_nor, b1, b2, b3, b4);

   % get values with upper change of a
   lce_nor(ia) = lce_nor(ia) - 2*delta_a;
   fact_dw = gaussianFunctionAct(lce_nor, b1, b2, b3, b4);

   lce_nor(ia) = lce_nor(ia) + delta_a;

   df_dlce_nor_fd(:, ia) = (fact_up - fact_dw)/(2*delta_a);

end

% check the drivative differences
tolerance = 1e-6;

% difference check of df_da
errorid_df_dlce_nor_fd = diffEvaluate(dfact_dlce_nor_equ, df_dlce_nor_fd, tolerance);

if ~isempty(errorid_df_dlce_nor_fd)
   sprintf('Differentiations in df_dlce_nor beyond thresholds\n')
end
