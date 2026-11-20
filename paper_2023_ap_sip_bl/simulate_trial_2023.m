function trial = simulate_trial_2023(params, SNR_db, trial_idx)
%SIMULATE_TRIAL_2023 Run one deterministic AP-SIP PA-BL trial.
%   Seeds rng(rng_seed_start + trial_idx - 1, 'threefry'). BER is NaN and detection
%   is skipped unless params.compute_ber (resolved by run_sim_2023) is true.
    local_validate_inputs(params, SNR_db, trial_idx);
    rng(params.rng_seed_start + trial_idx - 1, 'threefry');

    total_timer = tic;
    % SNR_db is the data-symbol SNR: n0 = data_power_2023 / 10^(SNR/10).
    n0 = params.data_power_2023 * 10^(-SNR_db / 10);
    chanParams = generate_channel_params(params);
    if strcmpi(params.paper_gain_mode, 'unity')
        chanParams.pathGains = ones(size(chanParams.pathGains));
    end
    truth = generate_paper_channel_2023(params, chanParams);
    frame = generate_ap_sip_frame(params);
    [Ydd, ~] = apply_sparse_otfs_channel_2023( ...
        frame.XDD, truth, params, n0);
    estimator_params = local_estimator_params(params);

    Yp = local_decouple(Ydd, frame.P);
    OmegaPilot = build_ap_sip_dictionary( ...
        frame.Xp, frame.P, params, 'pilot');
    pa_timer = tic;
    [Hpa, paState] = estimate_pa_bl_mimo(Yp, OmegaPilot, ...
        n0 * eye(size(OmegaPilot, 1)), estimator_params);
    paState.Hhat = Hpa;
    pa_runtime = toc(pa_timer);

    if params.compute_ber
        detection_timer = tic;
        [~, bitsHat] = detect_lmmse_2023(Ydd, Hpa, frame, estimator_params, n0);
        detection_runtime = toc(detection_timer);
        pa_ber = mean(bitsHat(:) ~= frame.bits(:));
    else
        detection_runtime = 0;
        pa_ber = NaN;
    end

    support_threshold = local_support_threshold(params);
    pa_metrics = compute_csi_metrics_2023( ...
        truth.H, Hpa, truth.support, paState.lambda, support_threshold, params);
    [true_parameters, estimated_parameters] = local_parameter_summary( ...
        chanParams, truth, Hpa, paState.lambda, params);

    trial = struct( ...
        'pa_nmse', pa_metrics.nmse, ...
        'pa_ber', pa_ber, ...
        'pa_support', pa_metrics, ...
        'pa_iterations', paState.iterations, ...
        'pa_converged', paState.converged, ...
        'pa_runtime_seconds', pa_runtime, ...
        'detection_runtime_seconds', detection_runtime, ...
        'true_parameters', true_parameters, ...
        'estimated_parameters', estimated_parameters, ...
        'runtime_seconds', toc(total_timer));
end

function observations = local_decouple(Ydd, decoder)
    observations = complex(zeros(size(Ydd, 1) * size(decoder, 2), size(Ydd, 3)));
    for rx = 1:size(Ydd, 3)
        observations(:, rx) = reshape(Ydd(:, :, rx) * decoder, [], 1);
    end
end

function threshold = local_support_threshold(params)
    threshold = 0.1;
    if isfield(params, 'support_threshold_2023')
        threshold = params.support_threshold_2023;
    end
    if ~isscalar(threshold) || ~isreal(threshold) || ~isfinite(threshold) || ...
            threshold < 0 || threshold > 1
        error('simulate_trial_2023:InvalidSupportThreshold', ...
            'support_threshold_2023 must be a finite real scalar in [0, 1].');
    end
end

function [true_parameters, estimated_parameters] = local_parameter_summary( ...
        chanParams, truth, Hhat, lambda, params)
    true_parameters = struct( ...
        'pathDelays', chanParams.pathDelays(:), ...
        'pathDopplers', chanParams.pathDopplers(:), ...
        'AoAs', chanParams.AoAs(:), ...
        'AoDs', chanParams.AoDs(:), ...
        'pathGains', truth.pathGains(:));

    Q = params.M_tau * params.G_nu;
    estimate_count = min(params.targets_num, numel(lambda));
    [~, ranked_atoms] = sort(lambda(:), 'descend');
    atom_indices = ranked_atoms(1:estimate_count);
    delay_indices = mod(atom_indices - 1, params.M_tau) + 1;
    doppler_indices = floor((atom_indices - 1) / params.M_tau) + 1;
    coefficient_norm = zeros(numel(atom_indices), 1);
    for index = 1:numel(atom_indices)
        rows = atom_indices(index) + (0:params.Nt - 1) * Q;
        coefficient_norm(index) = norm(Hhat(rows, :), 'fro');
    end
    estimated_parameters = struct( ...
        'atom_indices', atom_indices, ...
        'pathDelays', params.delay_grid_2023(delay_indices), ...
        'pathDopplers', params.doppler_grid_2023(doppler_indices), ...
        'lambda', lambda(atom_indices), ...
        'coefficient_norm', coefficient_norm);
end

function estimator_params = local_estimator_params(params)
    fields = {'M', 'N', 'Nt', 'Nr', 'K2', 'M_tau', 'G_nu', ...
        'delay_grid_2023', 'doppler_grid_2023', 'ModOrder', ...
        'em_max_iter_2023', 'em_tol_2023', ...
        'variance_floor_2023', 'data_power_2023'};
    estimator_params = struct();
    for index = 1:numel(fields)
        field = fields{index};
        if ~isfield(params, field)
            error('simulate_trial_2023:MissingEstimatorParameter', ...
                'params.%s is required by the AP-SIP estimators.', field);
        end
        estimator_params.(field) = params.(field);
    end
end


function local_validate_inputs(params, SNR_db, trial_idx)
    if ~isstruct(params)
        error('simulate_trial_2023:InvalidParams', 'params must be a struct.');
    end
    if ~isscalar(SNR_db) || ~isreal(SNR_db) || ~isfinite(SNR_db)
        error('simulate_trial_2023:InvalidSNR', 'SNR_db must be a finite real scalar.');
    end
    if ~isscalar(trial_idx) || ~isreal(trial_idx) || ~isfinite(trial_idx) || ...
            trial_idx < 0 || trial_idx ~= round(trial_idx)
        error('simulate_trial_2023:InvalidTrialIndex', ...
            'trial_idx must be a nonnegative integer scalar.');
    end
end
