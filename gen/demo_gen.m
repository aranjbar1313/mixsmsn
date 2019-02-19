arg(1).mu = [2,2];
arg(2).mu = [-2,2];
arg(3).mu = [-2,-2];
arg(1).Sigma = [1,0.5;0.5,1];
arg(2).Sigma = [1,0.5;0.5,1];
arg(3).Sigma = [1,0.5;0.5,1];
arg(1).shape = [1,1];
arg(2).shape = [-1,1];
arg(3).shape = [-1,-1];
arg(1).nu = 4;
arg(2).nu = 4;
arg(3).nu = 4;
y = rmmix(1e5, [0.2,0.3,0.5], "Skew_t", arg);
histogram2(y(:,1),y(:,2),[200,200],'FaceColor','flat')