%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is to generate muscle forces and joint moments based on the
% optimized muscle parameters.
%
% By: Huawei Wang
% Date: August 8, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mom, force] = muscleForceMoment_MPO(x, lmt, d, M, N, S, J, P)

    lce_opt = x(end - P + 1:end - P + M);
    lt_slack = x(end - P + M + 1:end - P + 2*M);
    theta0 = x(end - P + 2*M + 1:end - P + 3*M);
    Fmax = x(end - P + 3*M + 1:end);
    
    mom = zeros(sum(N), J);
    force = zeros(sum(N), M);
    
    for t = 1:length(N)
        sta_st = sum(N(1:t-1))*M*S;
        mea_st = sum(N(1:t-1));
        for n = 1:N(t)
            sta_st_n = (n-1)*M*S;
            
            lce = x(sta_st + sta_st_n + 2*M + 1: sta_st + sta_st_n + 3*M);
            lmt_mea = lmt(mea_st+n, :);
            d_mea = d(mea_st+n, :);
            
            % calculate muscle force and joint torques
            force(mea_st+n, :) = tendenForce_Groote_MPO(lmt_mea, lce, lce_opt, lt_slack, theta0);
            mom(mea_st+n, :) = force(mea_st+n, :).*Fmax*d_mea';

        end
    end
end