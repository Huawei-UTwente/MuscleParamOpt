%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is to test the activation dynamics of Friedl's paper
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% some default parameters, from Groote 2016 paper
    
    b = 0.1; % transition smoothness parameter
    Ta = 0.015;  % activation time constant
    Td = 0.060;  % deactivation time constant
    
    delta_t = 0.01; % time interval of forward simulation
    
    % muscle excitation consist of several step changes
    u = [zeros(1, 50) + 0.5, zeros(1, 50) + 1, zeros(1, 50) + 0.5, ...
        zeros(1, 50) + 0, zeros(1, 50) + 0.5, zeros(1, 50) + 0, ...
        zeros(1, 50) + 1, zeros(1, 50) + 0];
    
    % initilize a and da
    a = zeros(1, length(u));
    a(1) = u(1);
    da = zeros(1, length(u));   

    % forward simulation of muscle activation (forward euler)
    
    for i = 1:length(u) - 1
        
    % equality constraints of the neural activation
    ft = 0.5*tanh(b*(u(i) - a(i)));
    da(i) = (1./(Ta.*(0.5 + 1.5*a(i))).*(ft + 0.5) + (0.5 + 1.5*a(i))./Td.*(-ft + 0.5)).*(u(i) - a(i));
    a(i+1) = a(i) + da(i)*delta_t;
    
    end

    newcolors = [0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
         
    colororder(newcolors)
    
    figure;
    st = 1;
    ed = length(lce);
    
    subplot(1,2,1)
    plot((1:length(u))*delta_t, u, 'k-', 'linewidth', 2);
    hold on
    plot((1:length(a))*delta_t, a, 'r--', 'linewidth', 2);
    legend(["u", "a"])
    
    subplot(1,2,2)
    plot((1:length(da))*delta_t, da, 'k-', 'linewidth', 2);
    
    
    
    
