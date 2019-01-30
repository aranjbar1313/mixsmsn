%% Density/CDF of SNC with scale location %%
function dens = dmvSNC (y, mu, Sigma, lambda, nu)

    % y: must be a matrix where each row has 
    %    a multivariate vector of dimension 
    %    ncol(y) = p. nrow(y) = sample size
    % mu, lambda: must be of the vector type of 
    %             the same dimension equal to ncol(y) = p
    % lambda: 1 x p
    % Sigma: Matrix p x p
    nrow = @(x) size(x,1); 
    ncol = @(x) size(x,2);
    n = nrow(y);
    p = ncol(y);
    addpath('utils');
    dens = 2 .* (nu(1) .* transpose(mvnpdf(y, mu, Sigma./nu(2))) .* normcdf(sqrt(nu(2)) .* ...
                                                                            transpose(sum(transpose(reshape(repmat(lambda * inv(matrix_sqrt(Sigma)), 1, n), p, n)) .* ...
                                                                                          (y - transpose(reshape(repmat(mu', 1, n), p, n))), 2))) + ...
                 (1 - nu(1)) .* transpose(mvnpdf(y, mu, Sigma)) .* normcdf(transpose(sum(transpose(reshape(repmat(lambda * inv(matrix_sqrt(Sigma)), 1, n), p, n)) .* ...
                                                                                         (y - transpose(reshape(repmat(mu', 1, n), p, n))), 2))));
    rmpath('utils');
end