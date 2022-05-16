%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Muscle dynamics and it's derivatives, take CEINMS manual as reference
%
% By: Huawei Wang
% Date: August 1, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, df_da, df_dda, df_du, df_dC1, df_dC2] =...
    activation_dyn_anne(a, da, u, C1, C2, M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the activation dynamics from the CEINMS manual
%   
% da = (u - a)*(c1*u + c2)  % delay is not
% considered right now, will add later on.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % coefficient relationship
    
    f = da - (u - a).*(C1*u + C2);  % equality constraints of the neural activation
    
    df_da = (C1*u + C2);
    df_dda = ones(1, M);
    df_du = -(C1*u + C2) - (u - a)*C1;
    df_dC1 = -(u - a).*u;
    df_dC2 = -(u - a);
 
end