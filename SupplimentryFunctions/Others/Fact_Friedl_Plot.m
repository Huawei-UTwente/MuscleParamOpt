%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is to plot the muscle force-length relationship of Friedl's paper
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% muscle model curves, parameters from Groote 2016 paper

%% muscle force-length curves

    b11 = 0.815;
    b21 = 1.055;
    b31 = 0.162;
    b41 = 0.063;
    b12 = 0.433;
    b22 = 0.717;
    b32 = -0.030;
    b42 = 0.200;
    b13 = 0.100;
    b23 = 1.000;
    b33 = 0.354;
    b43 = 0.000;
    
    d1 = -0.318;
    d2 = -8.149;
    d3 = -0.374;
    d4 = 0.886;

    % normalize muscle fiber lengths
    lce_nor = [-2 -1 0 0.5 0.7 0.9 0.99 0.998 0.999 1 1.05 1.1 1.2 1.3 1.4 1.5 1.6 1.601 1.602 1.8 2];
    % [-5 -4 -3 -2 -1 0 0.1 0.2 0.3 0.4 0.402 0.4035 0.52725 0.55 0.6 0.62875 0.65 0.7 0.71875 0.75 0.8 0.86125 0.9 0.95 0.98 1.045 1.05 1.1 1.2175 1.3 1.43875 1.5 1.6 1.7 1.8 1.9 2.0 2.5 3 4 5];

    % active force length relationship and its differentiatons
    fce1 = gaussianFunctionAct(lce_nor, b11, b21, b31, b41);
    fce2 = gaussianFunctionAct(lce_nor, b12, b22, b32, b42);
    fce3 = gaussianFunctionAct(lce_nor, b13, b23, b33, b43);
    
    fce = fce1 + fce2 + fce3;

    %% passive force length relationship and its differentiations
    kpe = 4.0;
    e0 = 0.6;
    
    fpee1 = exp(kpe*(lce_nor - 1)/e0);

    fpee = (fpee1 - 1)./(exp(kpe) - 1);
    
    fpee_pos = (sqrt(fpee.^2 + 1e-5) + fpee)/2;
    
    %% muscle force-velocity curves
    
        % force velocity relationship and its differentiations.
    
    dlce_nor = [-10 -7 -5 -3 -1 -0.8 -0.6 -0.3 -0.1 0 0.1 0.3 0.6 0.8 1 3 5 7 10];

    fv_logfun = (d2*dlce_nor/10 + d3) + sqrt(((d2*dlce_nor/10 + d3).^2) + 1);
    fv_log = log(fv_logfun);

    fv = d1*fv_log + d4;
    
    %% Tendon force calculation
    
    c1 = 0.2;
    c2 = 1; %0.995;
    c3 = 0.2; %0.25;
    kT = 17.5;

    lt_nor = [ -0.05 -0.02 -0.01 -0.005 -0.002 -0.001 0 0.00131 0.00281 0.00431 0.00581 0.00731 ... 
               0.00881 0.0103 0.0118 0.0123 0.02 0.03 0.04 0.05 0.06 0.08 0.1 0.125 0.15 0.175 0.2 0.225 0.25 0.275 0.3 0.325 0.35 0.375 0.4];
    
    Fse_nor = c1*exp(kT.*(lt_nor + 1 - c2)) - c3;
    
    % lt_nor = 0.98:0.0001:1.12;
    % Fse_nor = c1*exp(kT.*(lt_nor - c2)) - c3;
    % Fse_nor_pos = (sqrt(Fse_nor.^2 + 1e-5) + Fse_nor)/2;

%     figure
%     plot(lt_nor, Fse_nor)
%     hold on
    
    %% plot all curves
    
    
    newcolors = [0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
         
    colororder(newcolors)
    
    figure;
    
    subplot(1,3,1)
    plot(lt_nor, Fse_nor, '-','Color', [.7 .7 .7], 'linewidth', 2);
    hold on
    plot(lt_nor, Fse_nor_pos, 'r--', 'linewidth', 2);
    hold off
    xlabel('Lt-nor')
    ylabel('Ft')
    
    subplot(1,3,2)
    st = 1;
    ed = length(lce_nor);
    plot(lce_nor(st:ed), fce(st:ed), 'k-', 'linewidth', 2);
    hold on;
    plot(lce_nor(st:ed), fpee(st:ed), '-', 'Color', [.7 .7 .7], 'linewidth', 2);
    hold on
    plot(lce_nor(st:ed), fpee_pos(st:ed), 'r--', 'linewidth', 2);
    hold on
    plot(lce_nor(st:ed), fce(st:ed) + fpee(st:ed), '-', 'linewidth', 2);
    hold on
    plot(lce_nor(st:ed), fce(st:ed) + fpee_pos(st:ed), '--', 'linewidth', 2);
    hold off
    xlim([lce_nor(st), lce_nor(ed)]);
    ylim([-0.05, 1.5]);
    legend(["Guassion 1", "Guassion 2", "Guassion 3", "F-lce", "F-pas", "F-tol"])
    xlabel('Lce-nor')
    ylabel('Fact, Fpas')
    
    subplot(1,3,3)
    plot(dlce_nor, fv, '-','Color', [.7 .7 .7], 'linewidth', 2);
    xlabel('dLce-nor')
    ylabel('Fv')