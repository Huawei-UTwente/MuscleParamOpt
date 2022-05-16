%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Muscle dynamics and it's derivatives, take CEINMS manual as reference
%
% By: Huawei Wang
% Date: August 1, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, df_da, df_dda, df_du] =...
    activationDyn_Groote_diff(a, da, u, M)
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
    
    % derivitives
    dft_da = -0.5*(1 - tanh(b*(u - a)).^2).*b;
    dft_du = 0.5*(1 - tanh(b*(u - a)).^2).*b;
    
    df_dda = ones(1, M);
    df_da = (1./(Ta.*(0.5 + 1.5*a)).*(ft + 0.5) + (0.5 + 1.5*a)./Td.*(-ft + 0.5)) ...
          - (-1.5./(Ta.*(0.5 + 1.5*a).^2).*(ft + 0.5) + 1./(Ta.*(0.5 + 1.5*a)).*dft_da ...
          + 1.5./Td.*(-ft + 0.5) + (0.5 + 1.5*a)./Td.*(-dft_da)).*(u - a);
    
    df_du = - (1./(Ta.*(0.5 + 1.5*a)).*(ft + 0.5) + (0.5 + 1.5*a)./Td.*(-ft + 0.5)) ...
          - (1./(Ta.*(0.5 + 1.5*a)).*dft_du + (0.5 + 1.5*a)./Td.*(-dft_du)).*(u - a);
    
%     df_dTa = (1./(Ta.^2.*(0.5 + 1.5*a)).*(ft + 0.5)).*(u - a);
%     df_dTd = ((0.5 + 1.5*a)./Td.^2.*(-ft + 0.5)).*(u - a);
 
end