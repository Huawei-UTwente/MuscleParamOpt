function [row, col] = jacobianstructure_ipopt_rc_MPO(auxdata)

    lmt = auxdata.lmt;
    N = auxdata.N;
    M = auxdata.M;
    S = auxdata.S;
    C = auxdata.C;
    P = auxdata.P;
    
    t_em = auxdata.t_em;
    hs = auxdata.hs;
    
    [row, col] = jacobianStructure_MPO(M, S, C, N, P, lmt, t_em, hs);

end