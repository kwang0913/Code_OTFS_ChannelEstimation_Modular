function [Hhat, XdHat, state] = estimate_da_bl_mimo(Yjoint, frame, paState, params, n0)
%ESTIMATE_DA_BL_MIMO Alternate AP-SIP channel and soft-data posteriors.
%   paState.Hhat must be populated by the caller with the PA-BL Hhat output
%   immediately after estimate_pa_bl_mimo.  The estimator alternates a joint
%   dictionary/posterior update with an LMMSE soft-data update and uses the
%   PA-BL row variances as its shared hyperparameter initialization.

    [Q, Nt, Nr, max_iterations, tolerance, variance_floor] = ...
        local_validate_inputs(Yjoint, frame, paState, params, n0);
    Ydd = local_reconstruct_observation(Yjoint, frame, params);
    Hhat = paState.Hhat;
    lambda = max(paState.lambda, variance_floor);
    [XdHat, bitsHat] = detect_lmmse_2023(Ydd, Hhat, frame, params, n0, []);
    regularization = max(n0, eps);
    force_lazy = isfield(params, 'force_lazy_da_2023') && params.force_lazy_da_2023;
    use_lazy = force_lazy || Q * Nt > 1024 || params.M * params.N * Q * Nt > 1e7;
    posterior_mode = 'dense';
    if use_lazy
        posterior_mode = 'lazy_exact';
    end
    converged = false;
    variance_change = inf;
    for iteration = 1:max_iterations
        [Hhat, Sigma] = local_da_posterior(Yjoint, XdHat, frame, lambda, Nt, ...
            regularization, variance_floor, params, use_lazy);
        lambda_new = local_update_lambda(Hhat, Sigma, Q, Nt, Nr, variance_floor);
        if any(~isfinite(lambda_new))
            error('estimate_da_bl_mimo:NonFiniteState', ...
                'DA-BL variance update produced non-finite values.');
        end
        variance_change = norm(lambda_new - lambda)^2;
        lambda = lambda_new;
        [XdHat, bitsHat] = detect_lmmse_2023(Ydd, Hhat, frame, params, n0, Sigma);
        if variance_change <= tolerance, converged = true; break; end
    end
    [Hhat, Sigma] = local_da_posterior(Yjoint, XdHat, frame, lambda, Nt, ...
        regularization, variance_floor, params, use_lazy);
    [XdHat, bitsHat] = detect_lmmse_2023(Ydd, Hhat, frame, params, n0, Sigma);
    if any(~isfinite(Hhat(:))) || any(~isfinite(XdHat(:))) || any(~isfinite(lambda))
        error('estimate_da_bl_mimo:NonFiniteState', 'DA-BL posterior contains non-finite values.');
    end
    state = struct('lambda', lambda, 'Sigma', Sigma, 'iterations', iteration, ...
        'converged', converged, 'variance_change', variance_change, 'bitsHat', bitsHat, ...
        'posterior_mode', posterior_mode);
end

function [Hhat, Sigma] = local_da_posterior(Yjoint, Xd, frame, lambda, Nt, n0, floor_value, params, use_lazy)
    symbols = struct('Xd', Xd, 'Xp', frame.Xp);
    precoder = struct('D', frame.D, 'P', frame.P);
    if use_lazy
        descriptor = local_joint_lazy_descriptor(symbols, precoder, params);
        [Hhat, Sigma] = local_lazy_posterior(descriptor, Yjoint, lambda, Nt, n0, floor_value, params);
    else
        Omega = build_ap_sip_dictionary(symbols, precoder, params, 'joint');
        [Hhat, Sigma] = local_posterior(Omega, Yjoint, lambda, Nt, n0, floor_value);
    end
end

function descriptor = local_joint_lazy_descriptor(symbols, precoder, params)
% Create the same joint operator as the builder without materializing it.
    XDD = complex(zeros(params.M, params.N, params.Nt));
    for tx = 1:params.Nt
        XDD(:, :, tx) = symbols.Xd(:, :, tx) * precoder.D' + ...
            symbols.Xp(:, :, tx) * precoder.P';
    end
    descriptor = struct('kind', 'ap_sip_lazy', 'mode', 'joint', 'XDD', XDD, ...
        'D', precoder.D, 'P', precoder.P, 'row_count', params.M * params.N, ...
        'column_count', params.M_tau * params.G_nu * params.Nt);
end

function Ydd = local_reconstruct_observation(Yjoint, frame, params)
    M = params.M;
    data_columns = params.N - params.K2;
    Ydd = complex(zeros(params.M, params.N, params.Nr));
    data_count = M * data_columns;
    for rx = 1:params.Nr
        Yd = reshape(Yjoint(1:data_count, rx), M, data_columns);
        Yp = reshape(Yjoint(data_count + 1:end, rx), M, params.K2);
        Ydd(:, :, rx) = Yd * frame.D' + Yp * frame.P';
    end
end

function [Hhat, Sigma] = local_posterior(Omega, Yjoint, lambda, Nt, n0, variance_floor)
    Q = numel(lambda);
    prior_precision = kron(eye(Nt), diag(1 ./ max(lambda, variance_floor)));
    posterior_precision = (Omega' * Omega) / n0 + prior_precision;
    posterior_precision = (posterior_precision + posterior_precision') / 2;
    [posterior_factor, posterior_status] = chol(posterior_precision, 'lower');
    if posterior_status ~= 0
        error('estimate_da_bl_mimo:PosteriorNotPD', ...
            'DA-BL posterior precision must be Hermitian positive definite.');
    end
    Sigma = posterior_factor' \ (posterior_factor \ eye(Q * Nt));
    Sigma = (Sigma + Sigma') / 2;
    Hhat = posterior_factor' \ (posterior_factor \ ((Omega' * Yjoint) / n0));
end

function [Hhat, Sigma] = local_lazy_posterior(descriptor, Yjoint, lambda, Nt, n0, variance_floor, params)
% Stream the analytic Gram through delay-pair correlations.  The joint
% decoder is unitary, so this is exactly A^H*A without an MN-by-QNt matrix.
    Q = numel(lambda);
    coefficient_count = Q * Nt;
    gram = local_analytic_joint_gram(descriptor, params);
    rhs = complex(zeros(coefficient_count, size(Yjoint, 2)));
    for receiver = 1:size(Yjoint, 2)
        rhs(:, receiver) = apply_ap_sip_dictionary_2023(descriptor, ...
            Yjoint(:, receiver), params, 'adjoint');
    end
    prior_precision = kron(eye(Nt), diag(1 ./ max(lambda, variance_floor)));
    posterior_precision = gram / n0 + prior_precision;
    posterior_precision = (posterior_precision + posterior_precision') / 2;
    [posterior_factor, posterior_status] = chol(posterior_precision, 'lower');
    if posterior_status ~= 0
        error('estimate_da_bl_mimo:PosteriorNotPD', ...
            'DA-BL lazy posterior precision must be Hermitian positive definite.');
    end
    Sigma = posterior_factor' \ (posterior_factor \ eye(coefficient_count));
    Sigma = (Sigma + Sigma') / 2;
    Hhat = posterior_factor' \ (posterior_factor \ (rhs / n0));
end

function gram = local_analytic_joint_gram(descriptor, params)
% Exact A^H*A from delay-pair correlations; rows are never jointly materialized.
    M = params.M;
    N = params.N;
    Nt = params.Nt;
    Mtau = params.M_tau;
    Gnu = params.G_nu;
    Q = Mtau * Gnu;
    samples = (0:M-1).';
    phases = exp(1j * 2 * pi * ...
        (reshape(samples, M, 1, 1) - reshape(params.delay_grid_2023(:), 1, Mtau, 1)) .* ...
        reshape(params.doppler_grid_2023(:), 1, 1, Gnu) / (M * N));
    shifts = complex(zeros(M, N, Nt, Mtau));
    for transmitter = 1:Nt
        for delay_index = 1:Mtau
            shifts(:, :, transmitter, delay_index) = circshift(descriptor.XDD(:, :, transmitter), ...
                params.delay_grid_2023(delay_index), 1);
        end
    end
    gram = complex(zeros(Q * Nt));
    for transmitter_i = 1:Nt
        for delay_i = 1:Mtau
            phase_i = reshape(phases(:, delay_i, :), M, Gnu);
            indices_i = (transmitter_i - 1) * Q + (0:Gnu-1) * Mtau + delay_i;
            shifted_i = shifts(:, :, transmitter_i, delay_i);
            for transmitter_j = 1:Nt
                for delay_j = 1:Mtau
                    phase_j = reshape(phases(:, delay_j, :), M, Gnu);
                    indices_j = (transmitter_j - 1) * Q + (0:Gnu-1) * Mtau + delay_j;
                    correlation = sum(conj(shifted_i) .* shifts(:, :, transmitter_j, delay_j), 2);
                    gram(indices_i, indices_j) = phase_i' * (correlation .* phase_j);
                end
            end
        end
    end
    gram = (gram + gram') / 2;
end

function lambda = local_update_lambda(Hhat, Sigma, Q, Nt, Nr, variance_floor)
    lambda = zeros(Q, 1);
    transmitter_offsets = (0:Nt-1) * Q;
    for q = 1:Q
        rows = q + transmitter_offsets;
        mean_energy = sum(abs(Hhat(rows, :)).^2, 'all') / (Nt * Nr);
        posterior_variance = real(sum(diag(Sigma(rows, rows)))) / Nt;
        lambda(q) = max(mean_energy + posterior_variance, variance_floor);
    end
end


function [Q, Nt, Nr, max_iterations, tolerance, variance_floor] = ...
        local_validate_inputs(Yjoint, frame, paState, params, n0)
    required_params = {'M', 'N', 'Nt', 'Nr', 'K2', 'M_tau', 'G_nu', ...
        'delay_grid_2023', 'doppler_grid_2023', 'ModOrder', 'da_max_iter_2023', ...
        'em_tol_2023', 'variance_floor_2023'};
    for i = 1:numel(required_params)
        if ~isfield(params, required_params{i})
            error('estimate_da_bl_mimo:MissingParameter', ...
                'params.%s is required.', required_params{i});
        end
    end
    dimensions = [params.M, params.N, params.Nt, params.Nr, params.K2, ...
        params.M_tau, params.G_nu];
    if any(~isfinite(dimensions)) || any(dimensions < 1) || ...
            any(dimensions ~= round(dimensions)) || params.K2 >= params.N
        error('estimate_da_bl_mimo:InvalidParams', ...
            'Dimensions must be positive integers and K2 must be smaller than N.');
    end
    Nt = params.Nt;
    Nr = params.Nr;
    Q = params.M_tau * params.G_nu;
    max_iterations = params.da_max_iter_2023;
    tolerance = params.em_tol_2023;
    variance_floor = params.variance_floor_2023;
    if ~isscalar(max_iterations) || ~isreal(max_iterations) || ...
            ~isfinite(max_iterations) || max_iterations < 1 || ...
            max_iterations ~= round(max_iterations) || ...
            ~isscalar(tolerance) || ~isreal(tolerance) || ~isfinite(tolerance) || tolerance < 0 || ...
            ~isscalar(variance_floor) || ~isreal(variance_floor) || ...
            ~isfinite(variance_floor) || variance_floor <= 0
        error('estimate_da_bl_mimo:InvalidParams', ...
            'DA-BL iteration and variance parameters are invalid.');
    end
    if ~isnumeric(Yjoint) || ~isequal(size(Yjoint), [params.M * params.N, Nr]) || ...
            any(~isfinite(Yjoint(:)))
        error('estimate_da_bl_mimo:InvalidObservation', ...
            'Yjoint must stack vec(Y*D) then vec(Y*P) for every receiver.');
    end
    if ~isstruct(frame) || ~isfield(frame, 'D') || ~isfield(frame, 'P') || ...
            ~isfield(frame, 'Xp') || ~isfield(frame, 'bits') || ...
            ~isnumeric(frame.D) || ~isequal(size(frame.D), [params.N, params.N - params.K2]) || ...
            ~isnumeric(frame.P) || ~isequal(size(frame.P), [params.N, params.K2]) || ...
            ~isnumeric(frame.Xp) || ~isequal(size(frame.Xp), [params.M, params.K2, Nt]) || ...
            any(~isfinite(frame.D(:))) || any(~isfinite(frame.P(:))) || any(~isfinite(frame.Xp(:)))
        error('estimate_da_bl_mimo:InvalidFrame', ...
            'frame must provide finite AP-SIP D/P precoders and pilot symbols.');
    end
    if ~isstruct(paState) || ~isfield(paState, 'lambda') || ~isfield(paState, 'Hhat') || ...
            ~isnumeric(paState.lambda) || ~isequal(size(paState.lambda), [Q, 1]) || ...
            ~isnumeric(paState.Hhat) || ~isequal(size(paState.Hhat), [Q * Nt, Nr]) || ...
            any(~isfinite(paState.lambda)) || any(paState.lambda < 0) || ...
            any(~isfinite(paState.Hhat(:)))
        error('estimate_da_bl_mimo:InvalidPAState', ...
            'paState must provide finite lambda and Hhat outputs from PA-BL.');
    end
    if ~isscalar(n0) || ~isreal(n0) || ~isfinite(n0) || n0 < 0
        error('estimate_da_bl_mimo:InvalidNoiseVariance', ...
            'n0 must be a finite nonnegative real scalar.');
    end
end
