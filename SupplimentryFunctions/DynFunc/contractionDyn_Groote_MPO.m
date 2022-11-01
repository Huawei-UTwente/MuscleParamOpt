
function f = contractionDyn_Groote_MPO(lmt, a, lce, dlce, lce_opt, lt_slack, theta0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the contraction dynamics of Hill's muscle model
%
% By: Huawei Wang
% Date: August 1, 2020
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
    lce_nor = lce./lce_opt;

    % active force length relationship and its differentiatons
    fce = gaussianFunctionAct(lce_nor, b11, b21, b31, b41) + ...
          gaussianFunctionAct(lce_nor, b12, b22, b32, b42) + ...
          gaussianFunctionAct(lce_nor, b13, b23, b33, b43);

    % passive force length relationship and its differentiations
    fpee1 = exp(kpe*(lce_nor - 1)/e0);

    fpee = (fpee1 - 1)./(exp(kpe) - 1);
    
    fpee_s = (sqrt(fpee.^2 + 1e-5) + fpee)/2;

    % force velocity relationship and its differentiations.
    dlceMax = 10*lce_opt;
    dlce_nor = dlce./dlceMax;

    fv_logfun = (d2*dlce_nor + d3) + sqrt(((d2*dlce_nor + d3).^2) + 1);
    fv_log = log(fv_logfun);

    fv = d1*fv_log + d4;

    cos_theta = pennationAngSmooth(lce, lce_opt, theta0);
    
    % calcualte the force of the contraction element and PEE together

    Fce = (a.*fce.*fv + fpee_s).*cos_theta;

    % tendon force calculation
    Fse = tendenForce_Groote_MPO(lmt, lce, lce_opt, lt_slack, theta0);

    % force of the contraction element should equal to the force of the tense unit
    f = (Fce - Fse);

end 