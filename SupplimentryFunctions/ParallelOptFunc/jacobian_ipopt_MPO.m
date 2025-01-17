
function J = jacobian_ipopt_MPO(x, auxdata)

    lmt = auxdata.lmt;
    N = auxdata.N;
    M = auxdata.M;
    S = auxdata.S;
    C = auxdata.C;
    P = auxdata.P;
    t_em = auxdata.t_em;
    hs = auxdata.hs;
    
    row = auxdata.row;
    col = auxdata.col;

    
    jac_nz = jacobian_MPO(x, M, S, C, N, P, lmt, t_em, hs);

    J = sparse(row, col, jac_nz, auxdata.rJac, auxdata.cJac);

end