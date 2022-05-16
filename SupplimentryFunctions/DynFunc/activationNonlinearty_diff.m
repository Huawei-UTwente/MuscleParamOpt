function dadu = activationNonlinearty_diff(u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function determines the nonlinearty of the muscle activation, from
% the recorded EMG signals. To be simplified, a linear relationship can be
% used, however, to achieve more percise muscle forces and joint torques,
% the nonlinear mapping method from Kurt Manal can be used:
% https://www.sciencedirect.com/science/article/pii/S0021929003001520?via%3Dihub

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% if a simple linear mapping is neede, then:

% a = u;

%%
% otherwise, a nonlinear mapping is defined as below;

A = 0.1;  % the main parameter that determine the curvture, can be adjusted
u0 = 0.3085 - A*cos(pi/4);  % transition point;
a0 = 0.3085 + A*sin(pi/4); % transition point;

% determine the parameters for the curvture, according to the paper
m = (1 - a0)/(1-u0);
c = 1 - m;

%% sloving the alpha and beta parameters through an iteration process. Only
% need to solve once and then using the solved values (below) for the
% optimization

% alpha = solvingAlpha(a0, u0, m);
% beta = (exp(a0/alpha) - 1)/u0;
% 
% % check the nonlinearty curvture, based on the evaluated parameters
% u = 0:0.01:1;
% a = zeros(size(u));
% for i = 1:length(u)
%    if u(i) < u0
%        a(i) = alpha*log(beta*u(i) + 1);
%    else
%        a(i) = m*u(i) + c;
%    end
% end
% plot(u, a, '-*')
% axis equal
% grid
% hold off

%%
alpha = 0.2467;
beta = 15.3524;

dadu(u < u0) = alpha./(beta*u(u < u0) + 1).*beta;
dadu(u >= u0) = m;

end