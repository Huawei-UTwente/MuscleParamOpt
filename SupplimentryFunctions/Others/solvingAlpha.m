function alpha_i = solvingAlpha(a0, u0, m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to solve the aplha value of the activation nonlearty, 
% given the continuious condition.
% accroding to: https://www.sciencedirect.com/science/article/pii/S0021929003001520?via%3Dihub
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha_i_1 = 1;
for i = 1:10000
    [fa, dfa] = alphaEqualityConstraint(alpha_i_1, a0, u0, m);
    alpha_i = alpha_i_1 - fa/dfa;
    while abs(alpha_i - alpha_i_1) < 1e-4
        break
    end
    alpha_i_1 = alpha_i;
end
end


function [fa, dfa] = alphaEqualityConstraint(x, a0, u0, m)
    beta = (exp(a0/x) - 1)/u0;
    fa = m - x*beta/(beta*u0 + 1);
    dfa = 1/u0*(-1 + (1- a0/x)*exp(-a0/x));
end