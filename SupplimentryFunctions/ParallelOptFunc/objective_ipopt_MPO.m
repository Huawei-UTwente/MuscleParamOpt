function f = objective_ipopt_MPO(x, auxdata)

    lmt = auxdata.lmt;
    torque = auxdata.torque;
    mus_act = auxdata.mus_act;
    d = auxdata.d;
    muscle_par0 = auxdata.muscle_par0;
    N = auxdata.N;
    M = auxdata.M;
    S = auxdata.S;
    P = auxdata.P;
    W1 = auxdata.W1;
    W2 = auxdata.W2;
    W3 = auxdata.W3;
    W4 = auxdata.W4;
    W5 = auxdata.W5;

    % calculate objective function
    f =  objective_MPO(x, lmt, torque, mus_act, d, M, N, S, P,...
                          muscle_par0, W1, W2, W3, W4, W5);

end