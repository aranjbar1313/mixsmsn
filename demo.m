% demo
clear;clc;
addpath(genpath('.'))
arg = [-3,2,0,10;3,3,0,10];
[y,clu] = rmix(1e6, [0.3,0.7], 't', arg);
% histogram(y,1000)
settings.family = 't';
settings.get_init = true;
settings.criteria = true;
settings.group = false;
settings.error = 1e-5;
settings.iter_max = 100;
settings.calc_im = true;
settings.obs_prob = false;
initial_values.g = 2;
init_nu = 9;

out = smsn_mix(y, init_nu , initial_values, settings);
