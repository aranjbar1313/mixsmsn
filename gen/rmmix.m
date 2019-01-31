function [x1,clu] = rmmix(n, pii, family, arg)
% Function to generate mixtures of g multivariate distributions
% n: number of samples generated
% p: mixtures proportion vector (size g)
% family: Family of distribution
% arg: must be an struct array. Fields are arguments pass to the generator
% function. Fields include: n, mu, Sigma, shape, nu

if((family ~= "t") && (family ~= "Skew_t") && (family ~= "Skew_cn") && (family ~= "Skew_slash") && (family ~= "Skew_normal") && (family ~= "Normal")) 
    error("Family "+family + " Not Recognized. (List of Families (case-sensitive)): t, Skew_t , Skew_cn, Skew_slash, Skew_normal, Normal");
end


if (family == "Normal")
    rF1 = "gen_SN_multi";
    for i = 1 : length(arg)
        arg(i).shape = zeros(1,length(arg(i).mu));
        if (length(fieldnames(arg(i))) ~= 4 && length(fieldnames(arg(i))) ~= 3)
            error("Some arguments are missing for number of argument = " + num2str(i))
        end
    end
end


if (family == "Skew_normal")
    rF1 = "gen_SN_multi";
    for i = 1:length(arg)
        if (length(fieldnames(arg(i))) ~= 4 && length(fieldnames(arg(i))) ~= 3)
            error("Some arguments are missing for number of argument = " + num2str(i))
        end
    end
end


if ((family == "t") || (family == "Skew_t"))
    rF1 = "gen_ST_multi";
    for i = 1 : length(arg)
        if (length(fieldnames(arg(i))) ~= 4)
            error("Some arguments are missing for number of argument = " + num2str(i))
        end
    end
    if (family == "t")
        for i = 1 : length(arg)
            arg(i).shape = zeros(1,length(arg(i).mu));
        end
    end
end


if (family == "Skew_cn")
    rF1 = "gen_SCN_multi";
    for i = 1 : length(arg)
        if (length(fieldnames(arg(i))) ~= 4)
            error("Some arguments are missing for number of argument = " + num2str(i))
        end
    end
end


if (family == "Skew_slash")
    rF1 = gen_SS_multi;
    for i = 1 : length(arg)
        if (length(fieldnames(arg(i))) ~= 4)
            error("Some arguments are missing for number of argument = " + num2str(i))
        end
    end
end


x1 = zeros(n,length(arg(1).mu));
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
    arg(i).n = sum(clu == i);
    x1(clu == i,:) = feval(rF1,arg(i));
end