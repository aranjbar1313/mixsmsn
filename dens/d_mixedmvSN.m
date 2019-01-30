function dens = d_mixedmvSN (y, pi1, mu, Sigma, lambda)

    % y: the data matrix
    % pi1: must be of the vector type of dimension g
    % mu: must be of type list with g entries. Each entry in the list must be a vector of dimension p
    % Sigma: must be of type list with g entries. Each entry in the list must be an matrix p x p
    % lambda: must be of type list with g entries. Each entry in the list must be a vector of dimension p

    g = numel(pi1);
    dens = 0;
    for (j = 1 : g) 
        dens = dens + pi1(j) .* dmvSN(y, mu{j}, Sigma{j}, lambda{j});
    end
end