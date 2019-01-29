%% Density of SNCI Mixtures %%
function dens = d_mixedSNC (x, pi1, mu, sigma2, shape, nu)

    % x: the data vector
    % other parameters must be of the row-vector type of dimension g (qtd of mixtures)
    g = numel(pi1);
    dens = 0;
    for j = 1 : g 
        dens = dens + pi1(j) .* dSNC(x, mu(j), sigma2(j), shape(j), nu);
    end
end