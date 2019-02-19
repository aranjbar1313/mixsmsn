function x1 =  gen_SCN_multi(n, mu, Sigma, shape, nu)

if isstruct(n)
    mu = n.mu;
    Sigma = n.Sigma;
    shape = n.shape;
    nu = n.nu;
    n = n.n;
end

x1 = zeros(n,length(mu));
% for i = 1 : n
%     u = rand(1);
%     if (u < nu(1))
%         x1(i,:) = gen.SN.multi(1, mu, Sigma/nu(2), shape);
%     else
%         x1(i,:) = gen.SN.multi(1, mu, Sigma, shape);
%     end
% end

u = rand(1,n);
x1(u <= nu(1),:) = gen_SN_multi(sum(u < nu(1)), mu, Sigma/nu(2), shape);
x1(u > nu(1),:)= gen_SN_multi(sum(u > nu(1)), mu, Sigma, shape);