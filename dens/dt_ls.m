%% ST density/CDF with scale location %%
function dens = dt_ls (x, varargin)
    
    varargin = validate_input(varargin, [0, 1, 1, 4]);
    [loc, sigma2, shape, nu] = varargin{:};
    
    d = (x - loc)/sqrt(sigma2);
    dens = 2 .* tpdf(d, nu) ...
             .* tcdf(sqrt((1 + nu)./(d.^2 + nu)) .* d .* shape, 1 + nu) ...
             ./ sqrt(sigma2);
end
