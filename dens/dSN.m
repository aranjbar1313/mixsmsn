%% SN density/CDF with scale location %%
function dens = dSN (y, varargin)
    
    varargin = validate_input(varargin, [0, 1, 1]);
    [mu, sigma2, shape] = varargin{:};

    dens = 2 * normpdf(y, mu, sqrt(sigma2)) ...
             * normcdf(shape*(y - mu)/sqrt(sigma2));
end