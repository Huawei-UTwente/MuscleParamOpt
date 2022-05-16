function cos_theta = pennationAngSmooth(lce, lce_opt, theta0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the pennation angle cos_theta and corresponding derivative
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % smooth function parameter that determine the smoothness of the
    % changes
    smooth_delta = 1e-8;

    % angle of the contraction element
    h = lce_opt.*sin(theta0);
    
    % calcualte the ce length in the parallel direction
    y = lce.^2 - h.^2; 

    dy = (sqrt(y.^2 + smooth_delta) + y)/2;
    
    cos_theta = sqrt(dy)./lce;  % calcualte cos(theta) at current lce length

end