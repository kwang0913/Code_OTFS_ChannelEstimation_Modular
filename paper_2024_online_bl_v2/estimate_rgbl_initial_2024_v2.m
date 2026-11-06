function [Hhat, state] = estimate_rgbl_initial_2024_v2(R, Psi, Rw, params)
%ESTIMATE_RGBL_INITIAL_2024_V2 Initial row-group Bayesian channel estimate.
    [Q, Nt, Nr, max_iterations, tolerance, floor_value] = ...
        local_validate(R, Psi, Rw, params);
    lambda = ones(Q, 1);
    converged = false;
    squared_changes = zeros(max_iterations, 1);
    for iteration = 1:max_iterations
        posterior = local_posterior(R, Psi, Rw, lambda, Nt, floor_value);
        lambda_new = local_update_lambda( ...
            posterior.mean, posterior.diagonal, Q, Nt, Nr, floor_value);
        squared_changes(iteration) = norm(lambda_new-lambda)^2;
        lambda = lambda_new;
        if squared_changes(iteration) <= tolerance
            converged = true;
            break;
        end
    end
    posterior = local_posterior(R, Psi, Rw, lambda, Nt, floor_value);
    Hhat = posterior.mean;
    if any(~isfinite(Hhat(:))) || any(~isfinite(lambda)) || ...
            any(~isfinite(posterior.factor(:)))
        error('estimate_rgbl_initial_2024_v2:NonFiniteState', ...
            'Initial RGBL produced a non-finite posterior state.');
    end
    state = struct('Hhat', Hhat, 'lambda', lambda, ...
        'prior_diagonal', posterior.prior_diagonal, ...
        'covariance_diagonal', posterior.diagonal, ...
        'covariance_factor', posterior.factor, ...
        'iterations', iteration, 'converged', converged, ...
        'squared_changes', squared_changes(1:iteration));
end

function posterior = local_posterior(R, Psi, Rw, lambda, Nt, floor_value)
    prior_diagonal = repmat(max(real(lambda(:)), floor_value), Nt, 1);
    weighted_dictionary = Psi .* prior_diagonal.';
    output_covariance = weighted_dictionary*Psi' + Rw;
    output_covariance = (output_covariance+output_covariance')/2;
    [factor, status] = chol(output_covariance, 'lower');
    if status ~= 0
        error('estimate_rgbl_initial_2024_v2:PosteriorNotPD', ...
            'The output covariance is not positive definite.');
    end
    solved_observation = factor'\(factor\R);
    solved_dictionary = factor'\(factor\weighted_dictionary);
    covariance_factor = weighted_dictionary'/factor';
    diagonal = max(real(prior_diagonal - ...
        sum(conj(weighted_dictionary).*solved_dictionary, 1).'), 0);
    posterior = struct('mean', weighted_dictionary'*solved_observation, ...
        'diagonal', diagonal, 'prior_diagonal', prior_diagonal, ...
        'factor', covariance_factor);
end

function lambda = local_update_lambda(mean_value, diagonal, Q, Nt, Nr, floor_value)
    lambda = zeros(Q, 1);
    for atom = 1:Q
        rows = atom+(0:Nt-1)*Q;
        lambda(atom) = max( ...
            sum(abs(mean_value(rows, :)).^2, 'all')/(Nt*Nr) + ...
            sum(diagonal(rows))/Nt, floor_value);
    end
end

function [Q, Nt, Nr, max_iterations, tolerance, floor_value] = ...
        local_validate(R, Psi, Rw, params)
    if ~isnumeric(R) || ~ismatrix(R) || isempty(R) || any(~isfinite(R(:)))
        error('estimate_rgbl_initial_2024_v2:InvalidObservation', ...
            'R must be a finite nonempty matrix.');
    end
    if ~isnumeric(Psi) || ~ismatrix(Psi) || size(Psi, 1) ~= size(R, 1) || ...
            any(~isfinite(Psi(:)))
        error('estimate_rgbl_initial_2024_v2:InvalidDictionary', ...
            'Psi must be finite and have one row per observation.');
    end
    if ~isnumeric(Rw) || ~isequal(size(Rw), [size(R,1), size(R,1)]) || ...
            any(~isfinite(Rw(:))) || ...
            norm(Rw-Rw', 'fro') > 100*eps(max(1, norm(Rw, 'fro')))
        error('estimate_rgbl_initial_2024_v2:InvalidNoise', ...
            'Rw must be a finite Hermitian covariance matching R.');
    end
    [~, noise_status] = chol((Rw+Rw')/2, 'lower');
    if noise_status ~= 0
        error('estimate_rgbl_initial_2024_v2:InvalidNoise', ...
            'Rw must be positive definite.');
    end
    required = {'Nt','Nr','em_max_iter_2024','em_tol_2024', ...
        'variance_floor_2024'};
    if ~isstruct(params) || ~isscalar(params) || ~all(isfield(params, required))
        error('estimate_rgbl_initial_2024_v2:InvalidParams', ...
            'Missing initial RGBL parameters.');
    end
    Nt = params.Nt; Nr = params.Nr;
    if size(R, 2) ~= Nr || mod(size(Psi, 2), Nt) ~= 0
        error('estimate_rgbl_initial_2024_v2:InvalidDimensions', ...
            'Nt/Nr do not match Psi/R.');
    end
    Q = size(Psi, 2)/Nt;
    max_iterations = params.em_max_iter_2024;
    tolerance = params.em_tol_2024;
    floor_value = params.variance_floor_2024;
    if any(~isfinite([Nt,Nr,max_iterations,tolerance,floor_value])) || ...
            any([Nt,Nr,max_iterations] < 1) || ...
            any([Nt,Nr,max_iterations] ~= round([Nt,Nr,max_iterations])) || ...
            tolerance < 0 || floor_value <= 0
        error('estimate_rgbl_initial_2024_v2:InvalidParams', ...
            'Invalid RGBL dimensions or stopping parameters.');
    end
end