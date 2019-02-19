function out = smsn_mmix(y, nu, initial_values, settings)
    
    nrow = @(x) size(x,1); 
    ncol = @(x) size(x,2);
    
    if (~ismatrix(y))
        error('The response is not in a matrix farmat.\n')
    end
    if (ncol(y) < 2)
        error('For the univariate cases use the smsn_mix function.\n')
    end
    if (isfield(settings, 'family'))
        family = settings.family;
        if (~strcmp(family, 'Skew.t') && ~strcmp(family, 't') && ...
            ~strcmp(family, 'Skew.normal') && ~strcmp(family, 'Normal'))
            error(strcat('Family ', family, ' not recognized.\n'))
        end
    else
        settings.family = 'Skew.normal';
    end

    if (~isfield(initial_values, 'g') && ~isfield(initial_values, 'params'))
        error('The model is not specified correctly.\n')
    end
    if (isfield(settings, 'get_init'))
        if (settings.get_init == false)
            g = length(initial_values.params);
            pii = initial_values.pii;
            for i = 1 : length(initial_values.params) - 1
                if (~isfield(initial_values.params(i), 'mu') || ~isfield(initial_values.params(i), 'Sigma') || ...
                    ~isfield(initial_values.params(i), 'shape'))
                    error('Please provide mu, Sigma, and shape for all groups.')
                end
                if (size(initial_values.params(i).mu) ~= size(initial_values.params(i + 1).mu) || ...
                    size(initial_values.params(i).Sigma) ~= size(initial_values.params(i + 1).Sigma) || ...
                    size(initial_values.params(i).shape) ~= size(initial_values.params(i + 1).shape))
                    error('The size of the initial values are not compatibles.\n')
                end
            end

            if ((length(initial_values.params) ~= length(pii)))
                error('The size of the initial values are not compatibles.\n')
            end
            if (sum(pii) ~= 1)
                error('probability of pii does not sum to 1.\n')
            end
        end
    end
    if (exist('g', 'var'))
        if (g < 1)
            error('g must be greater than 0.\n')
        end
    end
    p = ncol(y);
    n = nrow(y);
    if (settings.get_init == true)
        key = true;
        if (~isfield(initial_values, 'g'))
            error('g is not specified correctly.\n')
        else
            g = initial_values.g;
        end
        if isfield(initial_values, 'params')
            key = false;
        end
        if ~key
            if length(initial_values.params) ~= g
                error('The size of the initial values are not compatibles.\n')
            end
        end
        

        k_iter_max = 50;
        n_start = 1;
        if (isfield(settings, 'kmeans_param'))
            kmeans_param = settings.kmeans_param;
            if (isfield(kmeans_param, 'iter_max')) 
                k_iter_max = kmeans_param.iter_max;
            end
            if (isfield(kmeans_param, 'n_start'))
                n_start = kmeans_param.n_start;
            end
        end
        if g > 1
            if key
                [init.idx, init.C, init.sumd, init.D] = kmeans(y, g, 'MaxIter', k_iter_max, 'Replicates', n_start);
            else
                mu = [];
                for i = 1 : length(initial_values.params)
                    mu = cat(1, mu, initial_values.params(i).mu);
                end
                [init.idx, init.C, init.sumd, init.D] = kmeans(y, mu, 'MaxIter', k_iter_max, 'Replicates', n_start);
            end

            init.size = hist(init.idx);
            init.size = init.size(init.size ~= 0);
            pii = init.size/numel(y);
            mu = cell(1, g); shape = cell(1, g); Sigma = cell(1, g);
            for j = 1 : g
                if(key)
                    mu{j} = init.C(j, :);
                    shape{j} = sign(sum((y(init.idx == j, :) - repmat(mu{j}, [nrow(y(init.idx == j, :)), 1])).^3, 1));
                    Sigma{j} = cov(y(init.idx == j, :));
                end
            end
        end
    end

    if strcmp(family, 't')

    elseif strcmp(family, 'Skew.t')
        delta = cell(1, g);
        Delta = cell(1, g);
        Gama = cell(1, g);

        for k = 1 : g
            delta{k} = shape{k} ./ sqrt(1 + shape{k} * transpose(shape{k}));
            Delta{k} = transpose(matrix_sqrt(Sigma{k}) * transpose(delta{k}));
            Gama{k} = Sigma{k} - transpose(Delta{k}) * Delta{k};
        end

        if settings.uni_Gama
            Gama_uni = plus(Gama{:}) / g;
            Gama(:) = {Gama_uni};
        end
        
        % mu_old = mu;
        Delta_old = Delta;
        % Gama_old = Gama;

        criterio = 1;
        count = 0;
        lkante = 1;
        while ((criterio > settings.error) && (count <= settings.iter_max))
            tic
            count = count + 1;
            tal = zeros(n, g);
            S1 = zeros(n, g);
            S2 = zeros(n, g);
            S3 = zeros(n, g);
            for j = 1 : g 
                Dif = y - repmat(mu{j}, [n, 1]);
                Mtij2 = 1./(1 + Delta{j} * (Gama{j} \ transpose(Delta{j})));
                Mtij = sqrt(Mtij2);
                mtuij = sum(repmat(Mtij2 .* (Delta{j} / Gama{j}), [n, 1]) .* Dif, 2);
                A = mtuij ./ Mtij;
                mahalanobis = @(x, mu, Sigma) diag((x - mu) / Sigma * transpose(x - mu));
                dj = mahalanobis(y, mu{j}, Sigma{j});

                E = (2 .* (nu).^(nu./2) .* gamma((p + nu + 1)./2) .* ((dj + nu + A.^2)).^(-(p + nu + 1)./2))./(gamma(nu./2) .* (sqrt(pi)).^(p + 1) .* sqrt(det(Sigma{j})) .* dmvt_ls(y, mu{j}, Sigma{j}, shape{j}, nu));
                u = ((4 .* (nu).^(nu./2) .* gamma((p + nu + 2)./2) .* (dj + nu).^(-(p + nu + 2)./2))./(gamma(nu./2) .* sqrt(pi.^p) .* sqrt(det(Sigma{j})) .* dmvt_ls(y, mu{j}, Sigma{j}, shape{j}, nu))) .* tcdf(sqrt((p + nu + 2) ./ (dj + nu)) .* A, p + nu + 2);

                d1 = dmvt_ls(y, mu{j}, Sigma{j}, shape{j}, nu);
                if sum(d1 == 0)
                    d1(d1 == 0) = 1/intmax;
                end
                d2 = d_mixedmvST(y, pii, mu, Sigma, shape, nu);
                if sum(d2 == 0)
                    d2(d2 == 0) = 1/intmax;
                end

                tal(:, j) = d1 .* pii(j)./d2;
                S1(:, j) = tal(:, j) .* u;
                S2(:, j) = tal(:, j) .* (mtuij .* u + Mtij .* E);
                S3(:, j) = tal(:, j) .* (mtuij.^2 .* u + Mtij2 + Mtij .* mtuij .* E);

                pii(j) = (1./n) .* sum(tal(:, j));

                mu{j} = sum(S1(:, j) .* y - S2(:, j) .* repmat(Delta_old{j}, [n, 1]), 1)./sum(S1(:, j));
                Dif = y - mu{j};
                Delta{j} = sum(S2(:, j) .* Dif, 1)./sum(S3(:, j));

                sum2 = zeros(p);
                for i = 1 : n
                    sum2 = sum2 + (S1(i, j) .* (transpose(y(i, :) - mu{j})) * (y(i, :) - mu{j}) - ...
                                   S2(i, j) .* (transpose(Delta{j}) * (y(i, :) - mu{j})) - ...
                                   S2(i, j) .* (transpose(y(i, :) - mu{j}) * (Delta{j})) + ...
                                   S3(i, j) .* (transpose(Delta{j}) * (Delta{j})));
                end

                Gama{j} = sum2 ./ sum(tal(:, j));

                if ~settings.uni_Gama
                    Sigma{j} = Gama{j} + transpose(Delta{j}) * Delta{j};
                    shape{j} = (Delta{j} / matrix_sqrt(Sigma{j})) / (1 - Delta{j} / Sigma{j} * transpose(Delta{j})).^(1/2);
                end
            end
            %{
            if settings.uni_Gama
                GS = 0;
                for j = 1 : g 
                    GS = GS %+ (tal(:, j) .* Gama{j});
                end
                Gama_uni = transpose(sum(transpose(GS), 2))
                for j = 1 : g 
                    Gama{j} = Gama_uni;
                    Sigma{j} = Gama{j} + transpose(Delta{j}) * Delta{j};
                    shape{j} = zeros(1, p);
                end
            end
            %}
            logvero_ST = @(nu) -1*sum(log(d_mixedmvST(y, pii, mu, Sigma, shape, nu)));
            options = optimset('TolX', 0.000001);
            nu = fminbnd(logvero_ST, 0, 100, options);
            pii(g) = 1 - (sum(pii) - pii(g));

            zero_pos = pii == 0;
            pii(zero_pos) = 1e-10;
            pii(pii == max(pii)) = max(pii) - sum(pii(zero_pos));

            lk = sum(log(d_mixedmvST(y, pii, mu, Sigma, shape, nu)));
            criterio = abs((lk./lkante) - 1);

            lkante = lk;
            % mu_old = mu;
            Delta_old = Delta;
            % Gama_old = Gama;
            toc
        end
        if settings.criteria
            [~, cl] = max(tal, [], 2);
            icl = 0;
            for j = 1 : g
                icl = icl + sum(log(pii(j) .* dmvt_ls(y(cl == j), mu{j}, Sigma{j}, shape{j}, nu)));
            end
        end
    elseif strcmp(family, 'Normal')

    elseif strcmp(family, 'Skew.Normal')

    end 
    out.params = repmat(struct('mu', 0, 'Sigma', 0, 'shape', 0), [1, g]);
    for j = 1 : g 
        out.params(j).mu = mu{j};
        out.params(j).Sigma = Sigma{j};
        out.params(j).shape = shape{j};
    end
    out.nu = nu;
    out.pii = pii;
end
