function y = gen_Skew_cn(n, mu, sigma2, shape, nu1, nu2)
% Function to Generate Contaminated Skew Normal
% n: number of values to be generated
% mu, sigma2, shape: location, scale and asymmetry respectively

% rmix(n, nu[1], gen.Skew.normal, list(c(mu,sigma2/nu[2],shape), c(mu,sigma2,shape)))

if ~isscalar(n)
    t = num2cell(n);
    [n, mu, sigma2, shape, nu1, nu2] = deal(t{:});
end

y = rmix(n, [nu1,1-nu1], "Skew_normal", [[mu,sigma2/nu2,shape]; [mu,sigma2,shape]]);