
function c = constraints_ipopt_MPO(x, auxdata)

    lmt = auxdata.lmt;
    M = auxdata.M;
    S = auxdata.S;
    C = auxdata.C;
    N = auxdata.N;
    P = auxdata.P;
    hs = auxdata.hs;
    
    % equality constraints
    c = constraints_MPO(x, M, S, C, N, P, lmt, t_em, hs);

end