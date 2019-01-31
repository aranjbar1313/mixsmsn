function d = mahalanobis(x,center,S)

d = (x-center)*inv(S)*transpose(x-center);