function [x1,clu] = rmix(n, pii, family, arg)
% Function to generate mixtures of g distributions
% n: number of samples generated
% p: mixtures proportion vector (size g)
% family: Family of distribution
% arg: must be a matrix with each row containing a vector of size of
% arguments to be passed to family distribution generator

if ((family ~= "t") && (family ~= "Skew_t") && (family ~= "Skew_cn") && (family ~= "Skew_slash") && (family ~= "Skew_normal") && (family ~= "Normal"))
    error("Family "+family + " Not Recognized. (List of Families (case-sensitive)): t, Skew_t , Skew_cn, Skew_slash, Skew_normal, Normal");
end

if((family == "Normal") || (family == "Skew_normal") )
    rF1 = "gen_Skew_normal";
    for i = 1 : size(arg,1)
        if(length(arg(i,:)) ~= 4 && length(arg(i,:)) ~= 3)
            error("Number of arguments is not correct for argument  " + num2str(i))
        end
    end
    if(family == "Normal")
        for i = 1 : size(arg,1)
            arg(i,3) = 0;
        end
    end
end


if ((family == "t") || (family == "Skew_t"))
    rF1 = "gen_Skew_t";
    for i = 1 : size(arg,1)
        if(length(arg(i,:)) ~= 4)
            error("Number of arguments is not correct for argument  " + num2str(i))
        end
    end
    if(family == "t")
        for i = 1 : size(arg,1)
            arg(i,3) = 0;
        end
    end
end


if (family == "Skew_cn")
    rF1 = "gen_Skew_cn";
    for i = 1 : size(arg,1)
        if (length(arg(i,:)) ~= 5)
            error("Number of arguments is not correct for argument  " + num2str(i))
        end
    end
end


if (family == "Skew_slash")
    rF1 = "gen_Skew_slash";
    for i = 1 : size(arg,1)
        if (length(arg(i,:)) ~= 4)
            error("Number of arguments is not correct for argument  " + num2str(i))
        end
    end
end


x1 = zeros(1,n);
% clu = zeros(1,n);
g = length(pii);
% interval = 0;
% for j = 1 : g-1
%     interval = [interval, interval(j) + pii(j)];
% end
% interval = [interval, 1]
interval = [0,cumsum(pii)];
interval(end) = 1;
% for i = 1:n
%     u = rand(1)
%     clu(i) = find(interval>=u,1);
%     x1(i) = feval(rF1,[1,arg(clu(i),:)]);
% end
cluster_proportion = diff(interval);
clu = randsample(g,n,true,cluster_proportion);
for i = 1 : g
    x1(clu == i) = feval(rF1,[sum(clu==i),arg(i,:)]);
end





