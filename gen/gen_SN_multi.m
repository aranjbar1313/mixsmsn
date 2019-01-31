function y = gen_SN_multi(n, mu, Sigma, shape)

if isstruct(n)
    mu = n.mu;
    Sigma = n.Sigma;
    shape = n.shape;
    n = n.n;
end

p = length(mu);
delta = (shape / (sqrt(1 + shape*transpose(shape))))';
mu = mu';
% y = zeros(n,p);
addpath('../utils')
% for i = 1 : n
%     y(i,:) = mu + matrix_sqrt(Sigma)*(delta*abs(randn(1)) + matrix_sqrt(eye(p) - delta*transpose(delta))*mvnrnd(zeros(1,p),eye(p),1)');
% end
y = repmat(mu,1,n) + matrix_sqrt(Sigma)*(repmat(delta,1,n).*repmat(abs(randn(1,n)),p,1)  ...
    + matrix_sqrt(eye(p) - delta*transpose(delta))*mvnrnd(zeros(1,p),eye(p),n)');
rmpath('../utils')
y = y';
