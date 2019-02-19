clear;clc;
addpath(genpath('.'))

arg(1).mu = [ 2,  2];
arg(2).mu = [-2,  2];
arg(3).mu = [-2, -2];
arg(1).Sigma = [1, 0.5; 0.5, 1];
arg(2).Sigma = [1, 0.5; 0.5, 1];
arg(3).Sigma = [1, 0.5; 0.5, 1];
arg(1).shape = [ 1,  1];
arg(2).shape = [-1,  1];
arg(3).shape = [-1, -1];
arg(1).nu = 4;
arg(2).nu = 4;
arg(3).nu = 4;
pii = [0.2, 0.3, 0.5];
y = rmmix(1e4, pii, "Skew_t", arg);

initial_values.g  = 3;
init_nu           = 1;
settings.family   = 'Skew.t';
settings.get_init = true;
settings.criteria = true;
settings.group    = false;
settings.error    = 1e-4;
settings.iter_max = 100;
settings.uni_Gama = false;
settings.calc_im  = true;
settings.obs_prob = true;
out = smsn_mmix(y, init_nu, initial_values, settings);

s = 100;
a = reshape(linspace(-5, 5, s), [s, 1]);
[x, y] = meshgrid(a);
mu1 = {arg(1).mu, arg(2).mu, arg(3).mu};
Sigma1 = {arg(1).Sigma, arg(2).Sigma, arg(3).Sigma};
shape1 = {arg(1).shape, arg(2).shape, arg(3).shape};
nu1 = arg(1).nu;
pii1 = pii;
dens = d_mixedmvST ([x(:), y(:)], pii1, mu1, Sigma1, shape1, nu1);
dens = reshape(dens, [s, s]);
figure
contour(x, y, dens);
title('original pdf contours')

mu2 = {out.params(1).mu, out.params(2).mu, out.params(3).mu};
Sigma2 = {out.params(1).Sigma, out.params(2).Sigma, out.params(3).Sigma};
shape2 = {out.params(1).shape, out.params(2).shape, out.params(3).shape};
nu2 = out.nu;
pii2 = out.pii;
dens = d_mixedmvST ([x(:), y(:)], pii2, mu2, Sigma2, shape2, nu2);
dens = reshape(dens, [s, s]);
figure
contour(x, y, dens);
title('estimated pdf contours')