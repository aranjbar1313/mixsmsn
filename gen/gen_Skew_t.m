function y = gen_Skew_t(n, mu, sigma2, shape, nu)
% Function to generate Skew-t
% n: number of values to be generated
% mu, sigma2, shape: location, scale and asymmetry respectively

if ~isscalar(n)
    t = num2cell(n);
    [n, mu, sigma2, shape, nu] = deal(t{:});
end

y = mu + (random('gamma',nu/2, 2/nu,1,n)).^(-1/2).*gen_Skew_normal(n, 0, sigma2, shape);

end

