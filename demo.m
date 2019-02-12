% demo
clear;clc;
addpath(genpath('.'))
% arg = [-2, 12, -1, 23; 1, 17, 8, 23];
arg = [-3, 4, -3; 3, 2, 3];
[y,clu] = rmix(1e6, [0.3,0.7], 'Normal', arg);
x = linspace(min(y), max(y), 1000);
% dens = d_mixedST(x, [0.3,0.7], arg(:,1), arg(:,2), arg(:,3), arg(1,4));
dens = d_mixedSN(x, [0.3,0.7], arg(:,1), arg(:,2), [0, 0]);
histogram(y,1000,'Normalization','pdf')
hold on
plot(x,dens,'LineWidth',2)
settings.family = 'Normal';
settings.get_init = true;
settings.criteria = true;
settings.group = false;
settings.error = 1e-6;
settings.iter_max = 1000;
settings.calc_im = true;
settings.obs_prob = false;
initial_values.g = 2;
init_nu = 1;

tic
out = smsn_mix(y, init_nu , initial_values, settings);
toc
dens_pred = d_mixedSN(x, out.pii, out.mu, out.sigma2, out.shape);
% dens_pred = d_mixedST(x, out.pii, out.mu, out.sigma2, out.shape, out.nu);
plot(x,dens_pred,'--','LineWidth',3)