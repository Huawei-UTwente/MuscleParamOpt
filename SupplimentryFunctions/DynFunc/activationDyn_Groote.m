%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Muscle dynamics and it's derivatives, take CEINMS manual as reference
%
% By: Huawei Wang
% Date: August 1, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f =  activationDyn_Groote(a, da, u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the activation dynamics from the Groote 2016 muscle model
% https://doi.org/10.1007/s10439-016-1591-9
% By: Huawei Wang
% Date: 12/06/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    b = 0.1; % transition smoothness parameter
    Ta = 0.015;
    Td = 0.060;
    
    % equality constraints of the neural activation
    ft = 0.5*tanh(b*(u - a));
    f = da - (1./(Ta.*(0.5 + 1.5*a)).*(ft + 0.5) + (0.5 + 1.5*a)./Td.*(-ft + 0.5)).*(u - a);
 
end