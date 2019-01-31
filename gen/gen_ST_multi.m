function y = gen_ST_multi(n, mu, Sigma, shape, nu)

if isstruct(n)
    mu = n.mu;
    Sigma = n.Sigma;
    shape = n.shape;
    nu = n.nu;
    n = n.n;
end

y = repmat(mu,n,1) ...
    + repmat(random('gamma',nu/2, 2/nu,n,1),1,length(mu)).^(-1/2).*...
    gen_SN_multi(n, zeros(1, length(mu)), Sigma, shape);
