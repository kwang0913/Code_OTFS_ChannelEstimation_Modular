function results = run_sim_2024_v2(params)
%RUN_SIM_2024_V2 Aggregate main-physical RGBL/ORGBL trials.
    local_validate_params(params);
    n_trials = params.trials;
    n_snr = numel(params.SNR_db_list);
    metrics = local_empty_metrics(n_trials, n_snr);
    trials = cell(n_trials, n_snr);
    use_parallel = isfield(params, 'use_parallel') && params.use_parallel && ...
        n_trials > 1 && local_parallel_available();
    for snr_index = 1:n_snr
        snr_db = params.SNR_db_list(snr_index);
        fprintf('--- 2024 v2 main-physical BL SNR = %g dB ---\n', snr_db);
        temporary = cell(n_trials, 1);
        if use_parallel
            parfor trial_index = 1:n_trials
                temporary{trial_index} = simulate_trial_2024_v2( ...
                    params, snr_db, trial_index);
            end
        else
            for trial_index = 1:n_trials
                temporary{trial_index} = simulate_trial_2024_v2( ...
                    params, snr_db, trial_index);
            end
        end
        for trial_index = 1:n_trials
            trial = temporary{trial_index};
            metrics = local_store(metrics, trial, trial_index, snr_index);
        end
        trials(:, snr_index) = temporary;
        fprintf('Main NMSE (RGBL -> ORGBL): %.2f -> %.2f dB\n', ...
            10*log10(mean(metrics.NMSE_RGBL(:, snr_index))), ...
            10*log10(mean(metrics.NMSE_ORGBL(:, snr_index))));
        if all(metrics.detection_skipped(:, snr_index))
            fprintf('BER (ORGBL) skipped for multi-trial run.\n\n');
        else
            fprintf('BER (ORGBL): %.4e\n\n', ...
                mean(metrics.BER_ORGBL(:, snr_index)));
        end
    end
    results = metrics;
    results.trials = trials;
    results.SNR_db_list = params.SNR_db_list;
    results.params = params;
end

function metrics = local_empty_metrics(n_trials, n_snr)
    metrics = struct( ...
        'NMSE_RGBL', zeros(n_trials, n_snr), ...
        'NMSE_ORGBL', zeros(n_trials, n_snr), ...
        'coefficient_NMSE_RGBL', zeros(n_trials, n_snr), ...
        'coefficient_NMSE_ORGBL', zeros(n_trials, n_snr), ...
        'BER_ORGBL', NaN(n_trials, n_snr), ...
        'support_precision_RGBL', zeros(n_trials, n_snr), ...
        'support_recall_RGBL', zeros(n_trials, n_snr), ...
        'rgbl_iterations', zeros(n_trials, n_snr), ...
        'rgbl_converged', false(n_trials, n_snr), ...
        'orgbl_converged', false(n_trials, n_snr), ...
        'orgbl_innovation_norm', zeros(n_trials, n_snr), ...
        'orgbl_innovation_logdet', zeros(n_trials, n_snr), ...
        'orgbl_gain_norm', zeros(n_trials, n_snr), ...
        'rgbl_runtime_seconds', zeros(n_trials, n_snr), ...
        'orgbl_runtime_seconds', zeros(n_trials, n_snr), ...
        'detection_runtime_seconds', zeros(n_trials, n_snr), ...
        'detection_skipped', false(n_trials, n_snr), ...
        'runtime_seconds', zeros(n_trials, n_snr));
end

function metrics = local_store(metrics, trial, row, column)
    metrics.NMSE_RGBL(row,column) = trial.NMSE_RGBL;
    metrics.NMSE_ORGBL(row,column) = trial.NMSE_ORGBL;
    metrics.coefficient_NMSE_RGBL(row,column) = trial.coefficient_NMSE_RGBL;
    metrics.coefficient_NMSE_ORGBL(row,column) = trial.coefficient_NMSE_ORGBL;
    metrics.BER_ORGBL(row,column) = trial.BER_ORGBL;
    metrics.support_precision_RGBL(row,column) = trial.support_RGBL.precision;
    metrics.support_recall_RGBL(row,column) = trial.support_RGBL.recall;
    metrics.rgbl_iterations(row,column) = trial.rgbl_iterations;
    metrics.rgbl_converged(row,column) = trial.rgbl_converged;
    metrics.orgbl_converged(row,column) = all(trial.orgbl_converged);
    metrics.orgbl_innovation_norm(row,column) = mean(trial.orgbl_innovation_norm);
    metrics.orgbl_innovation_logdet(row,column) = mean(trial.orgbl_innovation_logdet);
    metrics.orgbl_gain_norm(row,column) = mean(trial.orgbl_gain_norm);
    metrics.rgbl_runtime_seconds(row,column) = trial.rgbl_runtime_seconds;
    metrics.orgbl_runtime_seconds(row,column) = trial.orgbl_runtime_seconds;
    metrics.detection_runtime_seconds(row,column) = trial.detection_runtime_seconds;
    metrics.detection_skipped(row,column) = trial.detection_skipped;
    metrics.runtime_seconds(row,column) = trial.runtime_seconds;
end

function local_validate_params(params)
    if ~isstruct(params) || ~isscalar(params) || ...
            ~isfield(params, 'trials') || ~isfield(params, 'SNR_db_list') || ...
            ~isscalar(params.trials) || params.trials < 1 || ...
            params.trials ~= round(params.trials) || ...
            ~isnumeric(params.SNR_db_list) || isempty(params.SNR_db_list) || ...
            any(~isfinite(params.SNR_db_list(:)))
        error('run_sim_2024_v2:InvalidParams', ...
            'params must contain valid trials and SNR_db_list.');
    end
end

function available = local_parallel_available()
    available = false;
    try
        available = license('test', 'Parallel_Computing_Toolbox') || ...
            license('test', 'Distrib_Computing_Toolbox');
    catch
    end
end
