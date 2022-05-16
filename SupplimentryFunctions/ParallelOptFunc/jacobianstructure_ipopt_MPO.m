function J = jacobianstructure_ipopt_MPO(auxdata)

    row = auxdata.row;
    col = auxdata.col;
    
    N = auxdata.N;
    M = auxdata.M;
    S = auxdata.S;
    P = auxdata.P;

    J = sparse(row, col, ones(1, length(row)), sum(N-1)*C, ...
        sum(N)*M*S + P);

end