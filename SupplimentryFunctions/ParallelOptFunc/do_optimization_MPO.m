function do_optimization_MPO(opt, auxdata)

    %% generate initial guesses
    
    folder = auxdata.folder;
    options = auxdata.options;
    funcs = auxdata.funcs;
    
    M = auxdata.M;
    S = auxdata.S;
    N = auxdata.N;
    
    hs = auxdata.hs;
    mus_act = auxdata.mus_act;
    
    mus_par0 = auxdata.mus_par0;
    
    lb = options.lb;
    ub = options.ub;
    

    %% generate initial guesses
    % initialize the optimizing parameters
    mus_a = (0.95 + 0.1*rand(sum(N), M)).*mus_act(1:sum(N), :);
    mus_da = zeros(size(mus_a));
    for t = 1:length(N)
        if t== 1
            mus_da(2:N(t), :) = (mus_act(2:N(t), :) - mus_act(1:N(t)-1, :))./hs(t);
        else
            mus_da(sum(N(1:t-1))+2:sum(N(1:t)), :) = ...
                (mus_act(sum(N(1:t-1))+2:sum(N(1:t)), :) ...
                - mus_act(sum(N(1:t-1))+1:sum(N(1:t))-1, :))./hs(t);
        end
    end

    mus_da(1, :) = (0.95 + 0.1*rand(1, M)).*mus_da(2, :);
    mus_s = (0.9 + 0.2*rand(sum(N), M)).*mus_act(1:sum(N), :);

    lce = (0.95 + 0.1*rand(sum(N), M)).*mus_par0(1:M);
    dlce = (0.5 - 1*rand(sum(N), M)).*mus_par0(1:M);

    x0_1 = [mus_a(1:sum(N), :), mus_da(1:sum(N), :),...
        lce(1:sum(N), :), dlce(1:sum(N), :), mus_s];

    x0_2 = reshape(x0_1', [1, sum(N)*M*S]);
        
    par0 = lb(end - P + 1:end) ....
        + (ub(end - P + 1:end) - lb(end - P + 1:end)).*rand(1, P);

    x0 = [x0_2, par0];

    % The callback functions.
    funcs.objective         = @objective_ipopt_MPO;
    funcs.constraints       = @constraints_ipopt_MPO;
    funcs.gradient          = @gradient_ipopt_MPO;
    funcs.jacobian          = @jacobian_ipopt_MPO;
    funcs.jacobianstructure = @jacobianstructure_ipopt_MPO;

    [x, info] = ipopt(x0, funcs, options);

    saving_names = sprintf('%s/optimization_res%02d.mat', folder, opt);
    
    save_res_MPO(saving_names, x, info, auxdata)  % save optimized results
        
end