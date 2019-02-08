%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%             EM algorithm for mixing             %%%%%%%%%%%%%

function out = smsn_mix (y, nu, initial_values, settings)
    % y: the data vector (sample) of size n
    % mu, sigma2, shape, pii: are the initial values for the EM algorithm. 
    %   Each of them must be a vector of size g (the algorithm understands 
    %   the number of components to adjust based on the size of these vectors)
    % nu: initial value for the nu (in the case of the SNC must be 
    %   a two-dimensional vector with values between 0 and 1)
    % loglik: true of false - If you want to calculate the log-likelihood of the adjusted model
    % cluster: true or false - If you want a vector of size n indicating which component 
    %   the i-th observation belongs to
    % type: ['ST', 'SNC', 'SS', 'SN', 'Normal'] - ST: Adjusts by mixes of Skew-T
    %                                           - SNC: Adjusts by mixtures of Skew Normal Contaminated
    %                                           - SS: Adjusts for Skew-Slash Blends (ATENTION: STILL WITH PROBLEMS)
    %                                           - SN: Adjusts for Normal Skew Blends
    %                                           - Normal: Normal Mixes
    % error: Defines the algorithm stop criterion.

    nrow = @(x) size(x,1); 
    ncol = @(x) size(x,2);
%     mu = 0;

    if (nrow(y) > 1)
        error('This function is only for univariate response y!')
    end
    if (isfield(settings, 'family'))
        family = settings.family;
        if (strcmp(family, 't') && strcmp(family, 'Skew.t') && strcmp(family, 'Skew.cn') && ...
            strcmp(family, 'Skew.slash') && strcmp(family, 'Skew.normal') && strcmp(family, 'Normal'))
            error(strcat('Family ', family, ' not recognized.\n'))
        end
    else
        settings.family = 'Skew.normal';
    end
    if (~isfield(initial_values, 'g') && ((~isfield(initial_values, 'mu')) || (~isfield(initial_values, 'sigma2')) || ...
                                          (~isfield(initial_values, 'shape')) || (~isfield(initial_values, 'pii'))))
        error('The model is not specified correctly.\n')
    end
    if (isfield(settings, 'get_init'))
        if (settings.get_init == false)
            g = nrow(initial.mu);
            sigma2 = initial_values.sigma2;
            pii = initial_values.pii;
            if ((g ~= nrow(sigma2)) || (g ~= nrow(pii)))
                error('The size of the initial values are not compatibles.\n')
            end
        end
    else
        settings.get_init = true;
    end
    if (exist('g', 'var'))
        if (g < 1)
            error('g must be greater than 0.\n')
        end
    end
    
    if (settings.get_init == true)
        key.mu = true; key.sig = true; key.shp = true;
        if (~isfield(initial_values, 'g'))
            error('g is not specified correctly.\n')
        else
            g = initial_values.g;
            mu = zeros(1, g);
        end
        if (isfield(initial_values, 'mu'))
            if (nrow(initial_values.mu) > 0)
                key.mu = false;
            end
        end
        if (isfield(initial_values, 'sigma2'))
            if (nrow(initial_values.sigma2) > 0)
                key.sig = false;
            end
        end
        if (isfield(initial_values, 'shape'))
            if (nrow(initial_values.shape) > 0)
                key.shp = false;
            end
        end

        if (((~key.mu) && (nrow(initial_values.mu) ~= initial_values.g)) || ...
            ((~key.sig) && (nrow(initial_values.sigma2) ~= initial_values.g)) || ...
            ((~key.shp) && (nrow(initial_values.shape) ~= initial_values.g)))
            error('The size of the initial values are not compatibles.\n')
        end

        k_iter_max = 50;
        n_start = 1;
%         algorithm = 'Hartigan-Wong';
        if (isfield(settings, 'kmeans_param'))
            kmeans_param = settings.kmeans_param;
            if (isfield(kmeans_param, 'iter_max')) 
                k_iter_max = kmeans_param.iter_max;
            end
            if (isfield(kmeans_param, 'n_start'))
                n_start = kmeans_param.n_start;
            end
            if (isfield(kmeans_param, 'algorithm'))
%                 algorithm = kmeans_param.algorithm;
            end
        end
    
        if (g > 1)
            y = y';
            if (key.mu)
                [init.idx, init.C, init.sumd, init.D] = kmeans(y, g, 'MaxIter', k_iter_max, 'Replicates', n_start);
            else
                [init.idx, init.C, init.sumd, init.D] = kmeans(y, mu, 'MaxIter', k_iter_max, 'Replicates', n_start);
            end
            y = y';

            init.size = hist(init.idx);
            init.size = init.size(init.size ~= 0);
            pii = init.size/ncol(y);

            if (key.mu)
                mu = init.C;
                mu = mu';
            end
            if (key.sig)
                sigma2 = init.sumd' ./ init.size;
            end
            if (key.shp)
                shape = zeros(1, g);
                for j = 1 : g
                    m3 = (1 ./ init.size(j)) .* sum( (y(init.idx == j) - mu(j)).^3 );
                    shape(j) = sign(m3./(sigma2(j).^(3/2)));
                end
            end
        else
            if (key.mu)
                mu = mean(y);
            end
            if (key.sig)
                sigma2 = var(y);
            end
            pii = 1;
            if (key.shp)
                m3 = (1./numel(y)) * sum( (y - mu).^3 );
                shape = sign(m3./(sigma2.^(3/2)));
            end
        end
    end
    addpath('./dens')
    if strcmp(settings.family, 't')
        
        
        shape = zeros(1, g);
        lk = sum(log(d_mixedST(y, pii, mu, sigma2, shape, nu)));
        n = numel(y);
        Gama = zeros(1, g);
        Delta = Gama;
        delta = Gama;
        for (k = 1 : g)
            Gama(k) = sigma2(k) - Delta(k).^2;
        end
        teta = cat(2, mu, Delta, Gama, pii, nu);
        mu_old = mu;
        Delta_old = Delta;
        Gama_old = Gama;

        criterio = 1;
        count = 0;
        
        while ((criterio > settings.error) && (count <= settings.iter_max))
            
            count = count + 1;
            tal = zeros(g, n);
            S1 = zeros(g, n);
            S2 = zeros(g, n);
            S3 = zeros(g, n);
            for j = 1 : g
                
                tic
                dj = ((y - mu(j))./(sqrt(sigma2(j)))).^2;
                Mtij2 = 1./(1 + (Delta(j).^2)*(Gama(j).^(-1)));
                Mtij = sqrt(Mtij2);
                mutij = Mtij2 .* Delta(j) .* (Gama(j).^(-1)) .* (y - mu(j));
                A = mutij ./ Mtij;

                E = (2 .* (nu).^(nu./2) .* gamma((2 + nu)./2) .* ...
                    ((dj + nu + A.^2)).^(-(2 + nu)./2)) ./ (gamma(nu./2) .* pi .* sqrt(sigma2(j)) .* ...
                    dt_ls(y, mu(j), sigma2(j), shape(j) ,nu));
                u = ((4 .* (nu).^(nu./2) .* gamma((3 + nu)./2) .* (dj + nu).^(-(nu + 3)./2)) ./ ...
                    (gamma(nu./2) .* sqrt(pi) .* sqrt(sigma2(j)) .* ...
                    dt_ls(y, mu(j), sigma2(j),shape(j) ,nu)) ) .* ...
                    tcdf(sqrt((3 + nu)./(dj + nu)) .* A, 3+nu);
                
                d1 = dt_ls(y, mu(j), sigma2(j), shape(j), nu);
                if (numel(d1 == 0))
                    d1(d1 == 0) = 1/intmax;
                end
                d2 = d_mixedST(y, pii, mu, sigma2, shape, nu);
                if (numel(d2 == 0))
                    d2(d2 == 0) = 1/intmax;
                end

                tal(j, :) = d1 .* pii(j)./d2;
                S1(j, :) = tal(j, :) .* u;
                S2(j, :) = tal(j, :) .* (mutij .* u + Mtij .* E);
                S3(j, :) = tal(j, :) .* (mutij .^ 2 .* u + Mtij2 + Mtij .* mutij .* E);

                pii(j) = (1./n) .* sum(tal(j, :));
                mu(j) = sum(S1(j, :) .* y - Delta_old(j) .* S2(j, :))./sum(S1(j, :));
                Delta(j) = 0;
                Gama(j) = sum(S1(j, :) .* (y - mu(j)).^2 - 2 .* (y - mu(j)) .* Delta(j) .* S2(j, :) + Delta(j).^2 .* S3(j, :))./sum(tal(j, :));
                sigma2(j) = Gama(j) + Delta(j).^2;
                shape(j) = 0;
            end
            nbins = (max(y) - min(y)) * 100;
            [counts, centers] = hist(y, nbins);
            logvero_ST = @(nu) -1*sum(counts .* log(d_mixedST(centers, pii, mu, sigma2, shape, nu)));
            options = optimset('TolX', 0.000001);
            nu = fminbnd(logvero_ST, 0, 100, options);
            lk1 = sum(log(d_mixedST(y, pii, mu, sigma2, shape, nu)));
            pii(g) = 1 - (sum(pii) - pii(g));
            
            zero_pos = pii == 0;
            if (any(zero_pos))
                pii(zero_pos) = 10^-10;
                pii(pii == max(pii)) = max(pii) - sum(pii(zero_pos));
            end

            param = teta;
            teta = cat(2, mu, Delta, Gama, pii, nu);
            criterio = abs(lk1./lk - 1);

            mu_old = mu;
            Delta_old = Delta;
            Gama_old = Gama;
            lk = lk1;
            disp(count)
        end

        if (settings.criteria)
            [~, cl] = max(tal, [], 2);
            icl = 0;
            for j = 1 : g
                icl = icl + sum(log(pii(j) .* dt_ls(y(cl == j), mu(j), sigma2(j), shape(j), nu)));
            end
        end

        disp(count)
    elseif strcmp(settings.family, 'Skew.t')
        lk = sum(log(d_mixedST(y, pii, mu, sigma2, shape, nu)));
        n = numel(y);
        Gama = zeros(1, g);
        Delta = Gama;
        delta = Gama;
        for k = 1 : g
            delta(k) = shape(k) ./ (sqrt(1 + shape(k).^2));
            Delta(k) = sqrt(sigma2(k)) .* delta(k);
            Gama(k) = sigma2(k) - Delta(k).^2;
        end
        
        teta = cat(2, mu, Delta, Gama, pii, nu);
        mu_old = mu;
        Delta_old = Delta;
        Gama_old = Gama;

        criterio = 1;
        count = 0;

        while ((criterio > settings.error) && (count <= settings.iter_max))
            disp(count)
            count = count + 1;
            tal = zeros(g, n);
            S1 = zeros(g, n);
            S2 = zeros(g, n);
            S3 = zeros(g, n);
            for j = 1 : g
                dj = ((y - mu(j))./(sqrt(sigma2(j)))).^2;
                Mtij2 = 1./(1 + (Delta(j).^2)*(Gama(j).^(-1)));
                Mtij = sqrt(Mtij2);
                mutij = Mtij2 .* Delta(j) .* (Gama(j).^(-1)) .* (y - mu(j));
                A = mutij ./ Mtij;

                E = (2 .* (nu).^(nu./2) .* gamma((2 + nu)./2) .* ...
                    ((dj + nu + A.^2)).^(-(2 + nu)./2)) ./ (gamma(nu./2) .* pi .* sqrt(sigma2(j)) .* ...
                    dt_ls(y, mu(j), sigma2(j), shape(j) ,nu));
                u = ((4 .* (nu).^(nu./2) .* gamma((3 + nu)./2) .* (dj + nu).^(-(nu + 3)./2)) ./ ...
                    (gamma(nu./2) .* sqrt(pi) .* sqrt(sigma2(j)) .* ...
                    dt_ls(y, mu(j), sigma2(j),shape(j) ,nu)) ) .* ...
                    tcdf(sqrt((3 + nu)./(dj + nu)) .* A, 3+nu);

                d1 = dt_ls(y, mu(j), sigma2(j), shape(j), nu);
                if (numel(d1 == 0))
                    d1(d1 == 0) = 1/intmax;
                end
                d2 = d_mixedST(y, pii, mu, sigma2, shape, nu);
                if (numel(d2 == 0))
                    d2(d2 == 0) = 1/intmax;
                end

                tal(j, :) = d1 .* pii(j)./d2;
                S1(j, :) = tal(j, :) .* u;
                S2(j, :) = tal(j, :) .* (mutij .* u + Mtij .* E);
                S3(j, :) = tal(j, :) .* (mutij .^ 2 .* u + Mtij2 + Mtij .* mutij .* E);

                pii(j) = (1./n) .* sum(tal(j, :));
                mu(j) = sum(S1(j, :) .* y - Delta_old(j) .* S2(j, :))./sum(S1(j, :));
                Delta(j) = sum(S2(j, :) .* (y - mu(j))) ./ sum(S3(j, :));
                Gama(j) = sum(S1(j, :) .* (y - mu(j)).^2 - 2 .* (y - mu(j)) .* Delta(j) .* S2(j, :) + Delta(j).^2 .* S3(j, :))./sum(tal(j, :));
                sigma2(j) = Gama(j) + Delta(j).^2;
                shape(j) = ((sigma2(j).^(-1/2)) .* Delta(j))./(sqrt(1 - (Delta(j).^2).*(sigma2(j).^(-1))));
            end

            nbins = (max(y) - min(y)) * 100;
            [counts, centers] = hist(y, nbins);
            logvero_ST = @(nu) -1*sum(counts .* log(d_mixedST(centers, pii, mu, sigma2, shape, nu)));
%             logvero_ST = @(nu) -1*sum(log(d_mixedST(y, pii, mu, sigma2, shape, nu)));
            options = optimset('TolX', 0.000001);
            nu = fminbnd(logvero_ST, 0, 100, options);
            lk1 = sum(log(d_mixedST(y, pii, mu, sigma2, shape, nu)));
            pii(g) = 1 - (sum(pii) - pii(g));

            zero_pos = pii == 0;
            if (any(zero_pos))
                pii(zero_pos) = 10^-10;
                pii(pii == max(pii)) = max(pii) - sum(pii(zero_pos));
            end

            param = teta;
            teta = cat(2, mu, Delta, Gama, pii, nu);
            criterio = abs(lk1./lk - 1);

            mu_old = mu;
            Delta_old = Delta;
            Gama_old = Gama;
            lk = lk1;
        end

        if (settings.criteria)
            [~, cl] = max(tal, [], 2);
            icl = 0;
            for j = 1 : g
                icl = icl + sum(log(pii(j) .* dt_ls(y(cl == j), mu(j), sigma2(j), shape(j), nu)));
            end
        end
        

    elseif strcmp(family, 'Skew.normal')
        lk = sum(log(d_mixedSN(y, pii, mu, sigma2, shape)));
        n = numel(y);
        Gama = zeros(1, g);
        Delta = Gama;
        delta = Gama;
        for k = 1 : g 
            delta(k) = shape(k) ./ (sqrt(1 + shape(k).^2));
            Delta(k) = sqrt(sigma2(k)) .* delta(k);
            Gama(k) = sigma2(k) - Delta(k).^2;
        end

        teta = cat(2, mu, Delta, Gama, pii, nu);
        mu_old = mu;
        Delta_old = Delta;
        Gama_old = Gama;

        criterio = 1;
        count = 0;

        while ((criterio > settings.error) && (count <= settings.iter_max))
            count = count + 1;
            tal = zeros(g, n);
            S1 = zeros(g, n);
            S2 = zeros(g, n);
            S3 = zeros(g, n);
            for j = 1 : g
                Mtij2 = 1./(1 + (Delta(j).^2) .* (Gama(j).^(-1)));
                Mtij = sqrt(Mtij2);
                mutij = Mtij2 .* Delta(j) .* (Gama(j).^(-1)) .* (y - mu(j));

                prob = normcdf(mutij ./ Mtij);
                if (numel(prob == 0))
                    prob(prob == 0) = 1/intmax;
                end
                E = normpdf(mutij./mutij) ./ prob;
                u = ones(1, n);
                
                d1 = dSN(y, mu(j), sigma2(j), shape(j));
                if (numel(d1 == 0))
                    d1(d1 == 0) = 1/intmax;
                end
                d2 = d_mixedSN(y, pii, mu, sigma2, shape);
                if (numel(d2 == 0))
                    d2(d2 == 0) = 1/intmax;
                end

                tal(j, :) = d1 .* pii(j)./d2;

                S1(j, :) = tal(j, :) .* u;
                S2(j, :) = tal(j, :) .* (mutij .* u + Mtij .* E);
                S3(j, :) = tal(j, :) .* (mutij.^2 .* u + Mtij2 + Mtij .* mutij .* E);

                pii(j) = (1./n) .* sum(tal(j, :));
                mu(j) = sum(S1(j, :) .* y - Delta_old(j) .* S2(j, :))./sum(S1(j, :));
                Delta(j) = sum(S2(j, :) .* (y - mu(j)))./sum(S3(j, :));
                Gama(j) = sum(S1(j, :) .* (y - mu(j)).^2 - 2 .* (y - mu(j)) .* Delta(j) .* S2(j, :) + Delta(j).^2 * S3(j, :))./sum(tal(j, :));
                sigma2(j) = Gama(j) + Delta(j).^2;
                shape(j) = ((sigma2(j).^(-1/2)) .* Delta(j))./(sqrt(1 - (Delta(j).^2) .* (sigma2(j).^(-1))));
            end
            
            pii(g) = 1 - (sum(pii) - pii(g));

            zero_pos = pii == 0;
            if (numel(zero_pos))
                pii(zero_pos) = 10^-10;
                pii(pii == max(pii)) = max(pii) - sum(pii(zero_pos));
            end

            param = teta;
            teta = cat(2, mu, Delta, Gama, pii);
            lk1 = sum(log(d_mixedSN(y, pii, mu, sigma2, shape)));
            criterio = abs(lk1./lk - 1);

            mu_old = mu;
            Delta_old = Delta;
            Gama_old = Gama;
            lk = lk1;
        end

        if (settings.criteria)
            [~, cl] = max(tal, [], 2);
            icl = 0;
            for j = 1 : g
                icl = icl + sum(log(pii(j) .* dSN(y(cl == j), mu(j), sigma2(j), shape(j))));
            end
        end

    elseif strcmp(family, 'Normal')
        
        shape = zeros(1, g);
        lk = sum(log(d_mixedSN(y, pii, mu, sigma2, shape)));
        n = numel(y);
        Gama = zeros(1, g);
        Delta = Gama;
        delta = Gama;
        for k = 1 : g 
            Gama(k) = sigma2(k) - Delta(k).^2;
        end

        teta = cat(2, mu, Delta, Gama, pii);
        mu_old = mu;
        Delta_old = Delta;
        Gama_old = Gama;

        criterio = 1;
        count = 0;

        while ((criterio > settings.error) && (count <= settings.iter_max))
            count = count + 1;
            tal = zeros(g, n);
            S1 = zeros(g, n);
            S2 = zeros(g, n);
            S3 = zeros(g, n);
            for j = 1 : g
                Mtij2 = 1./(1 + (Delta(j).^2) .* (Gama(j).^(-1)));
                Mtij = sqrt(Mtij2);
                mutij = Mtij2 .* Delta(j) .* (Gama(j).^(-1)) .* (y - mu(j));

                prob = normcdf(mutij ./ Mtij);
                if (numel(prob == 0))
                    prob(prob == 0) = 1/intmax;
                end
                E = normpdf(mutij./mutij) ./ prob;
                u = ones(1, n);
                
                d1 = dSN(y, mu(j), sigma2(j), shape(j));
                if (numel(d1 == 0))
                    d1(d1 == 0) = 1/intmax;
                end
                d2 = d_mixedSN(y, pii, mu, sigma2, shape);
                if (numel(d2 == 0))
                    d2(d2 == 0) = 1/intmax;
                end

                tal(j, :) = d1 .* pii(j)./d2;
                S1(j, :) = tal(j, :) .* u;
                S2(j, :) = tal(j, :) .* (mutij .* u + Mtij .* E);
                S3(j, :) = tal(j, :) .* (mutij.^2 .* u + Mtij2 + Mtij .* mutij .* E);

                pii(j) = (1./n) .* sum(tal(j, :));
                mu(j) = sum(S1(j, :) .* y - Delta_old(j) .* S2(j, :))./sum(S1(j, :));
                Delta(j) = 0;
                Gama(j) = sum(S1(j, :) .* (y - mu(j)).^2 - 2 .* (y - mu(j)) .* Delta(j) .* S2(j, :) + Delta(j).^2 * S3(j, :))./sum(tal(j, :));
                sigma2(j) = Gama(j) + Delta(j).^2;
                shape(j) = 0;
            end
            
            pii(g) = 1 - (sum(pii) - pii(g));

            zero_pos = pii == 0;
            if (numel(zero_pos))
                pii(zero_pos) = 10^-10;
                pii(pii == max(pii)) = max(pii) - sum(pii(zero_pos));
            end

            param = teta;
            teta = cat(2, mu, Delta, Gama, pii);
            lk1 = sum(log(d_mixedSN(y, pii, mu, sigma2, shape)));
            criterio = abs(lk1./lk - 1);

            mu_old = mu;
            Delta_old = Delta;
            Gama_old = Gama;
            lk = lk1;
        end

        if (settings.criteria)
            [~, cl] = max(tal, [], 2);
            icl = 0;
            for j = 1 : g
                icl = icl + sum(log(pii(j) .* dSN(y(cl == j), mu(j), sigma2(j), shape(j))));
            end
        end

    end

    if settings.criteria
        if (strcmp(settings.family, 't') || strcmp(settings.family, 'Normal'))
            d = g .* 2 + (g - 1);
        else
            d = g .* 3 + (g - 1);
        end
        aic = -2 .* lk + 2 .* d;
        bic = -2 .* lk + log(n) .* d;
        edc = -2 .* lk + 0.2 .* sqrt(n) .* d;
        icl = -2 .* icl + log(n) .* d;
        
        out.aic = aic;
        out.bic = bic;
        out.edc = edc;
        out.icl = icl;
    end

    if ~settings.group

    end

    if settings.obs_prob
        
    end 

    if settings.calc_im

    end

    out.mu = mu;
    out.sigma2 = sigma2;
    out.shape = shape;
    out.nu = nu;
    out.pii = pii;
    out.iter = count;
    out.n = numel(y);
    out.group = cl;

end