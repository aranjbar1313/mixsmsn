function Asqrt = matrix_sqrt (A)

    [u, s, v] = svd(A);
    d = diag(s);
    if (min(d) >= 0)
        Asqrt = transpose(v * (transpose(u) .* sqrt(d)));
    else
        error('Matrix square root is not defined.\n')
    end
end