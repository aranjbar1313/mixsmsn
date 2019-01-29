%% Density of Skew Slash %%
function resp = dSS (y, mu, sigma2, shape, nu)

    resp = zeros(size(y));
    for i = 1 : numel(y)
        f = @(u) 2 .* nu .* u.^(nu - 1) ...
                   .* normpdf(y(i), mu, sqrt(sigma2./u)) ...
                   .* normcdf(u.^(1/2) .* shape .* (sigma2.^(-1/2)) .* (y(i) - mu));
        resp(i) = integral(f, 0, 1);
    end
end