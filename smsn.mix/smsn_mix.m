%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%             EM algorithm for mixing             %%%%%%%%%%%%%

function smsn_mix (y, nu, initial_values, settings)
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

    if (ncol(y) > 1)
        error('This function is only for univariate response y!')
    end
    if (isfield(settings, 'family'))
        family = settings.family;
        if ((family ~= 't') && (family ~= 'Skew.t') && (family ~= 'Skew.cn') && 
            (family ~= 'Skew.slash') && (family ~= 'Skew.normal') && (family ~= 'Normal'))
            error(strcat('Family ', family, ' not recognized.\n'))
        end
    else
        settings.family = 'Skew.normal';
    end
    if (!isfield(initial_values, 'g') && ((!isfield(initial_values, 'mu')) || (!isfield(initial_values, 'sigma2')) ||
                                          (!isfield(initial_values, 'shape')) || (!isfield(initial_values, 'pii'))))
        error('The model is not specified correctly.\n')
    end
    if (isfield(settings, 'get_init'))
        if (settings.get_init == false)
            g = nrow(initial.mu);
            if ((g ~= nrow(sigma2)) || (g ~= nrow(pii)))
                error('The size of the initial values are not compatibles.\n')
            end
        end
    else
        settings.get_init = true;
    end
    if (exist(g) && (g < 1))
        error('g must be greater than 0.\n')
    end
    if (settings.get_init == true)
        key.mu = true; key.sig = true; key.shp = true;
        if (!isfield(initial_values, 'g'))
            error('g is not specified correctly.\n')
        else
            g = initial_values.g;
        end
        if (nrow(initial_values.mu) > 0)
            key.mu = false;
        end
        if (nrow(initial_values.sigma2) > 0)
            key.sig = false;
        end
        if (nrow(initial_values.shape) > 0)
            key.shp = false;
        end

        if (((!key.mu) & (nrow(initial_values.mu) ~= initial_values.g)) |
            ((!key.sig) & (nrow(initial_values.sigma2) ~= initial_values.g)) |
            ((!key.shp) & (nrow(initial_values.shape) ~= initial_values.g)))
            error('The size of the initial values are not compatibles.\n')
        end

        k_iter_max = 50;
        n_start = 1;
        algorithm = 'Hartigan-Wong';
        if (isfield(settings, 'kmeans_param'))
            kmeans_param = settings.kmeans_param;
            if (isfield(kmeans_param, 'iter_max')) 
                k_iter_max = kmeans_param.iter_max;
            end
            if (isfield(kmeans_param, 'n_start'))
                n_start = kmeans_param.n_start;
            end
            if (isfield(kmeans_param, 'algorithm'))
                algorithm = kmeans_param.algorithm;
            end
        end
    end
    if (g > 1)
        if (key.mu)
            init = kmeans(y, g, 'MaxIter', k.iter_max, '')