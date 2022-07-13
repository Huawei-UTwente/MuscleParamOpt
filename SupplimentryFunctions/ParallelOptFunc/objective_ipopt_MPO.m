function f = objective_ipopt_MPO(x, auxdata)

    lmt = auxdata.lmt;
    torque = auxdata.torque;
    mus_act = auxdata.mus_act;
    d = auxdata.d;
    mus_par0 = auxdata.mus_par0;
    N = auxdata.N;
    M = auxdata.M;
    S = auxdata.S;
    P = auxdata.P;
    W1 = auxdata.W1;
    W2 = auxdata.W2;
    W3 = auxdata.W3;
    W4 = auxdata.W4;
    W5 = auxdata.W5;
    W6 = auxdata.W6;

    % calculate objective function
    f =  objective_MPO(x, lmt, torque, mus_act, d, M, N, S, P,...
                          mus_par0, W1, W2, W3, W4, W5, W6);

end