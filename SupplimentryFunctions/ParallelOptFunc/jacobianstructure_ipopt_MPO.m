function J = jacobianstructure_ipopt_MPO(auxdata)

    row = auxdata.row;
    col = auxdata.col;

    J = sparse(row, col, ones(1, length(row)), auxdata.rJac, auxdata.cJac);

end