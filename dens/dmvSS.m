%% Density/CDF of SS with scale location %%
function resp = dmvSS (y, mu, Sigma, lambda, nu)

    % y: must be a matrix where each row has 
    %    a multivariate vector of dimension 
    %    ncol(y) = p. nrow(y) = sample size
    % mu, lambda: must be of the vector type of 
    %             the same dimension equal to ncol(y) = p
    % lambda: 1 x p
    % Sigma: Matrix p x p
    mahalanobis = @(x, mu, Sigma) transpose(diag((x - mu) * inv(Sigma) * transpose(x - mu)));

    nrow = @(x) size(x,1); 
    ncol = @(x) size(x,2);
    n = nrow(y);
    p = ncol(y);
    resp = zeros(1, n);

    addpath('utils');
    for i = 1 : n
        di = mahalanobis(y(i,:), mu, Sigma);
        f = @(u) 2 .* nu .* u.^(nu - 1) .* ...
                    ((u./(2*pi)).^(p./2) .* det(Sigma).^(-1/2) .* exp(-u .* di./2)) .* ...
                    normcdf(u.^(1/2) .* (lambda * inv(matrix_sqrt(Sigma)) * transpose(y(i, :) - mu)));
        resp(i) = integral(f, 0, 1);
    end
    rmpath('utils');
end