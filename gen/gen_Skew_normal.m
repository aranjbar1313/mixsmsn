function y = gen_Skew_normal(n,mu,sigma2,shape)
% Function to generate random values of a Skew-Normal
% n: number of values to be generated
% mu, sigma2, shape: location, scale, and asymmetry respectively

if ~isscalar(n)
    t = num2cell(n);
    [n,mu,sigma2,shape] = deal(t{:});
end

delta = shape / sqrt(1 + shape^2);
y = mu*ones(1,n) + sqrt(sigma2)*(delta*abs(normrnd(0,1,1,n)) + ...
    (1 - delta^2)^(1/2)*normrnd(0,1,1,n));
end

