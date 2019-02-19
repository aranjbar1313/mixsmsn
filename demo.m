% demo
clear;clc;
addpath(genpath('.'))
% arg = [-2, 12, -1, 23; 1, 17, 8, 23];
arg = [-1, 4, -3; 1, 2, 3];
[y,clu] = rmix(1e5, [0.3,0.7], 'Skew_normal', arg);
x = linspace(min(y), max(y), 1000);
% dens = d_mixedST(x, [0.3,0.7], arg(:,1), arg(:,2), arg(:,3), arg(1,4));
dens = d_mixedSN(x, [0.3,0.7], arg(:,1), arg(:,2), arg(:, 3));
plot(x,dens,'LineWidth',2)
hold on
settings.family = 'Skew.normal';
settings.get_init = true;
settings.criteria = false;
settings.group = true;
settings.error = 1e-6;
settings.iter_max = 1000;
settings.calc_im = true;
settings.obs_prob = true;
initial_values.g = 2;
init_nu = 1;

tic
out = smsn_mix(y, init_nu , initial_values, settings);
toc
y1 = y(out.group == 1);
s1 = min(y1);
t1 = max(y1);
h = 1.1*max(dens);

rectangle('Position', [s1, 0, t1 - s1, h], 'FaceColor', [1, 0, 0, 0.2])
y2 = y(out.group == 2);
s2 = min(y2);
t2 = max(y2);
rectangle('Position', [s2, 0, t2 - s2, h], 'FaceColor', [0, 0, 1, 0.2])
histogram(y,300,'Normalization','pdf', 'FaceColor', 'k')
% histogram(y1, 100, 'Normalization', 'pdf', 'FaceColor', 'r')
% histogram(y2, 100, 'Normalization', 'pdf', 'FaceColor', 'b')
dens_pred = d_mixedSN(x, out.pii, out.mu, out.sigma2, out.shape);
% dens_pred = d_mixedST(x, out.pii, out.mu, out.sigma2, out.shape, out.nu);
plot(x,dens_pred,'--','LineWidth',3)