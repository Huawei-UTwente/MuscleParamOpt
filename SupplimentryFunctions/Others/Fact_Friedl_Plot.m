%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is to plot the muscle force-length relationship of Friedl's paper
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% some default parameters, from Groote 2016 paper
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
    kpe = 4.0;
    e0 = 0.6;
    d1 = -0.318;
    d2 = -8.149;
    d3 = -0.374;
    d4 = 0.886;

    % normalize muscle fiber lengths
    lce_nor = 0.0: 0.01: 2.0;

    % active force length relationship and its differentiatons
    fce1 = gaussianFunctionAct(lce_nor, b11, b21, b31, b41);
    fce2 = gaussianFunctionAct(lce_nor, b12, b22, b32, b42);
    fce3 = gaussianFunctionAct(lce_nor, b13, b23, b33, b43);
    
    fce = fce1 + fce2 + fce3;

    % passive force length relationship and its differentiations
    fpee1 = exp(kpe*(lce_nor - 1)/e0);

    fpee = (fpee1 - 1)./(exp(kpe) - 1);
    
    
    newcolors = [0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
         
    colororder(newcolors)
    
    figure;
    st = 40;
    ed = 160;
    plot(lce_nor(st:ed), fce1(st:ed), '--', 'linewidth', 2);
    hold on;
    plot(lce_nor(st:ed), fce2(st:ed), '--', 'linewidth', 2);
    hold on;
    plot(lce_nor(st:ed), fce3(st:ed), '--', 'linewidth', 2);
    hold on;
    plot(lce_nor(st:ed), fce(st:ed), 'k-', 'linewidth', 2);
    hold on;
    plot(lce_nor(st:ed), fpee(st:ed), '--', 'Color', [.7 .7 .7], 'linewidth', 2);
    xlim([lce_nor(st), lce_nor(ed)]);
    ylim([-0.05, 1.5]);
    legend(["Guassion 1", "Guassion 2", "Guassion 3", "F-lce", "F-pas"])
    xlabel('Lce-nor')
    ylabel('Fact, Fpas')
    
    
    