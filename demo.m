% demo
addpath(genpath('.'))
arg = [-3,2,0,17;3,3,0,10];
[x1,clu] = rmix(1e6, [0.3,0.7], "t", arg);
histogram(x1,1000)