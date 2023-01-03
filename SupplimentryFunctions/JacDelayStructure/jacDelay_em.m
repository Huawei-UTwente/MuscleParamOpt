%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function defined the effects of electromechancial delay on the jacobian
% structure
%
% By: Huawei Wang
% Date: April 14th, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function df_dx_cb = jacDelay_em(n_em, df_dx_em1, df_dx_em2, df_dx_em3, df_dx1,...
                             df_dx2, M, S)

    % depends on how large is the _em delay, x and x_em
    % maybe constructed using some common data nodes.
    % They should be combined together when extracting
    % the nonzero Jac elements.


    % if n_em is an integer multiple the time interval (hs)
    if ceil(n_em) == n_em  

        % beside the two df_dx, the combined fd_dx_cb length depends on the
        % n_em
        rblock = 2 + ceil(n_em);
        
        df_dx_cb = zeros(M, rblock*M*S);

        % df_dx are always at the end
        df_dx_cb(:, (rblock - 2)*M*S + 1:end) = [df_dx1, df_dx2];

        % df_dx_em1 does not exist, since stn_em1 = ceil(n_em)*M*S = floor(n_em)*M*S 
        df_dx_cb(:, 1:2*M*S) = df_dx_cb(:, 1:2*M*S) + [df_dx_em2, df_dx_em3];


    else  % if not, then the df_dx_em1 exist
        
        rblock = 3 + floor(n_em);

        df_dx_cb = zeros(M, rblock*M*S); % since three df_dx_em 

        % df_dx are always at the end
        df_dx_cb(:, (rblock - 2)*M*S + 1:end) = [df_dx1, df_dx2];

        % df_dx_em are always at the beginning, but now df_dx_em1 exist
        df_dx_cb(:, 1:3*M*S) = df_dx_cb(:, 1:3*M*S) + [df_dx_em1, df_dx_em2, df_dx_em3];

    end

end
