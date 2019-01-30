function ys = gen_Skew_slash(n, mu, sigma2, shape, nu)
% Function to generate Skew Slash
% n: number of values to be generated
% mu, sigma2, shape: location, scale and asymmetry respectively

if ~isscalar(n)
    t = num2cell(n);
    [n, mu, sigma2, shape, nu] = deal(t{:});
end

u1 = rand(1,n);
u2 = u1.^(1./nu);  % article formula (10) and inversion method
ys = mu + (u2).^(-1/2).*gen_Skew_normal(n, 0, sigma2, shape);