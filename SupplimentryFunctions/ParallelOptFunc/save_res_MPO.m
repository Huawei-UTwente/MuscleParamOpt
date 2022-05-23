function save_res_MPO(saving_names, x, info, auxdata)

    %% save results

    M = auxdata.M;
    S = auxdata.S;
    N = auxdata.N;
    J = auxdata.J;
    P = auxdata.P;
    lmt = auxdata.lmt;

    W1 = auxdata.W1;
    W2 = auxdata.W2;
    W3 = auxdata.W3;
    W4 = auxdata.W4;
    W5 = auxdata.W5;
    
    torque = auxdata.torque;
    mus_act = auxdata.mus_act;
    d = auxdata.d;
    mus_par0 = auxdata.mus_par0;

    % generate the joint moments and muscle force and activations

    [mom_res, force_res] = muscleForceMoment_MPO(x, lmt, d, M, N, S, J, P);

    obj = objective_MPO(x, lmt, torque, mus_act, d, M, N, S, P,...
                              mus_par0, W1, W2, W3, W4, W5);
    time = info.cpu;
    status = info.status;

    states = reshape(x(1:sum(N)*M*S), M*S, sum(N))';

    parameters = x(end - P + 1:end);

    save(saving_names, 'states', 'parameters', 'mom_res',...
        'force_res', 'obj', 'time', 'status');

end