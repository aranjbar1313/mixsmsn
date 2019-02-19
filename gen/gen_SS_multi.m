function ys =  gen_SS_multi(n, mu, Sigma, shape, nu)

if isstruct(n)
    mu = n.mu;
    Sigma = n.Sigma;
    shape = n.shape;
    nu = n.nu;
    n = n.n;
end

p = length(mu);
u1 = rand(n,1);
u2 = u1.^(1/(nu));   % article formula (10) and inversion method
ys = repmat(mu,n,1) + repmat((u2).^(-1/2),1,p).*gen_SN_multi(n, zeros(1,p), Sigma, shape);
