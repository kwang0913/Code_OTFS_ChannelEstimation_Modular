function results = run_sim_2023(params)
%RUN_SIM_2023 Aggregate deterministic AP-SIP PA-BL trials over SNR.
%   RESULTS = RUN_SIM_2023(PARAMS) returns PA-BL CSI and detection metrics
%   with rows for trials and columns for entries in PARAMS.SNR_db_list.

    local_validate_params(params);
    % Per-trial seeding switches the generator; restore the caller's state.
    caller_rng = rng;
    restore_rng = onCleanup(@() rng(caller_rng)); %#ok<NASGU>
    params.compute_ber = should_compute_ber(params);
    n_trials = params.trials;
    n_snr = numel(params.SNR_db_list);
    parallel_available = local_parallel_available();
    requested_parallel = isfield(params, 'use_parallel') && params.use_parallel;
    use_parallel = requested_parallel && n_trials > 1 && parallel_available;
    if requested_parallel && n_trials > 1 && ~parallel_available
        warning('run_sim_2023:ParallelUnavailable', ...
            'Parallel execution was requested but the Parallel Computing Toolbox is unavailable; using serial execution.');
    end

    NMSE_PA = zeros(n_trials, n_snr);
    BER_PA = zeros(n_trials, n_snr);
    support_precision = zeros(n_trials, n_snr);
    support_recall = zeros(n_trials, n_snr);
    pa_iterations = zeros(n_trials, n_snr);
    pa_converged = false(n_trials, n_snr);
    pa_runtime_seconds = zeros(n_trials, n_snr);
    detection_runtime_seconds = zeros(n_trials, n_snr);
    runtime_seconds = zeros(n_trials, n_snr);
    trials = cell(n_trials, n_snr);

    for snr_idx = 1:n_snr
        SNR_db = params.SNR_db_list(snr_idx);
        fprintf('--- 2023 AP-SIP PA-BL SNR = %g dB ---\n', SNR_db);

        temp_nmse = zeros(n_trials, 1);
        temp_ber = zeros(n_trials, 1);
        temp_precision = zeros(n_trials, 1);
        temp_recall = zeros(n_trials, 1);
        temp_iterations = zeros(n_trials, 1);
        temp_converged = false(n_trials, 1);
        temp_pa_runtime = zeros(n_trials, 1);
        temp_detection_runtime = zeros(n_trials, 1);
        temp_runtime = zeros(n_trials, 1);
        temp_trials = cell(n_trials, 1);
        local_update_progress('reset', n_trials);
        dq = local_progress_queue(parallel_available, n_trials);

        if use_parallel
            parallel_params = params;
            parallel_params.show_figures = false;
            parallel_params.verbose_output = false;
            parfor trial_idx = 1:n_trials
                trial = simulate_trial_2023(parallel_params, SNR_db, trial_idx);
                [temp_nmse(trial_idx), temp_ber(trial_idx), ...
                    temp_precision(trial_idx), temp_recall(trial_idx), ...
                    temp_iterations(trial_idx), temp_converged(trial_idx), ...
                    temp_pa_runtime(trial_idx), ...
                    temp_detection_runtime(trial_idx), temp_runtime(trial_idx)] = ...
                    local_trial_metrics(trial);
                temp_trials{trial_idx} = trial;
                local_send_progress(dq, n_trials);
            end
        else
            for trial_idx = 1:n_trials
                trial = simulate_trial_2023(params, SNR_db, trial_idx);
                [temp_nmse(trial_idx), temp_ber(trial_idx), ...
                    temp_precision(trial_idx), temp_recall(trial_idx), ...
                    temp_iterations(trial_idx), temp_converged(trial_idx), ...
                    temp_pa_runtime(trial_idx), ...
                    temp_detection_runtime(trial_idx), temp_runtime(trial_idx)] = ...
                    local_trial_metrics(trial);
                temp_trials{trial_idx} = trial;
                local_send_progress(dq, n_trials);
                if local_verbose(params)
                    fprintf('\n  Trial %d: PA-BL NMSE=%.2f dB, BER=%.4e\n', ...
                        trial_idx, 10 * log10(temp_nmse(trial_idx)), temp_ber(trial_idx));
                end
            end
        end
        fprintf('\n');

        NMSE_PA(:, snr_idx) = temp_nmse;
        BER_PA(:, snr_idx) = temp_ber;
        support_precision(:, snr_idx) = temp_precision;
        support_recall(:, snr_idx) = temp_recall;
        pa_iterations(:, snr_idx) = temp_iterations;
        pa_converged(:, snr_idx) = temp_converged;
        pa_runtime_seconds(:, snr_idx) = temp_pa_runtime;
        detection_runtime_seconds(:, snr_idx) = temp_detection_runtime;
        runtime_seconds(:, snr_idx) = temp_runtime;
        trials(:, snr_idx) = temp_trials;

        if n_trials == 1
            local_print_parameter_comparison(temp_trials{1});
        end
        if params.compute_ber
            fprintf('BER (PA-BL): %.4e\n', mean(temp_ber));
        else
            fprintf('BER (PA-BL) skipped (compute_ber = false).\n');
        end
        fprintf('NMSE (PA-BL operator): %.2f dB\n\n', ...
            10 * log10(mean(temp_nmse)));
    end

    results = struct( ...
        'NMSE_PA', NMSE_PA, ...
        'BER_PA', BER_PA, ...
        'support_precision', support_precision, ...
        'support_recall', support_recall, ...
        'pa_iterations', pa_iterations, ...
        'pa_converged', pa_converged, ...
        'pa_runtime_seconds', pa_runtime_seconds, ...
        'detection_runtime_seconds', detection_runtime_seconds, ...
        'runtime_seconds', runtime_seconds, ...
        'trials', {trials}, ...
        'SNR_db_list', params.SNR_db_list, ...
        'params', params);
end

function local_validate_params(params)
    if ~isstruct(params) || ~isfield(params, 'trials') || ~isfield(params, 'SNR_db_list')
        error('run_sim_2023:InvalidParams', ...
            'params must contain trials and SNR_db_list.');
    end
    if ~isscalar(params.trials) || ~isfinite(params.trials) || ...
            params.trials < 1 || params.trials ~= round(params.trials)
        error('run_sim_2023:InvalidTrials', 'params.trials must be a positive integer.');
    end
    if ~isvector(params.SNR_db_list) || isempty(params.SNR_db_list) || ...
            ~isreal(params.SNR_db_list) || any(~isfinite(params.SNR_db_list))
        error('run_sim_2023:InvalidSNRList', ...
            'params.SNR_db_list must be a nonempty finite real vector.');
    end
end

function available = local_parallel_available()
    available = false;
    try
        available = license('test', 'Parallel_Computing_Toolbox') || ...
            license('test', 'Distrib_Computing_Toolbox');
    catch
        available = false;
    end
end

function dq = local_progress_queue(parallel_available, total)
    dq = [];
    if parallel_available
        dq = parallel.pool.DataQueue;
        afterEach(dq, @(~) local_update_progress('increment', total));
    end
end

function local_send_progress(dq, total)
    if isempty(dq)
        local_update_progress('increment', total);
    else
        send(dq, 0);
    end
end

function [nmse, ber, precision, recall, iterations, converged, pa_runtime, ...
        detection_runtime, total_runtime] = local_trial_metrics(trial)
    nmse = trial.pa_nmse;
    ber = trial.pa_ber;
    precision = trial.pa_support.support_precision;
    recall = trial.pa_support.support_recall;
    iterations = trial.pa_iterations;
    converged = trial.pa_converged;
    pa_runtime = trial.pa_runtime_seconds;
    detection_runtime = trial.detection_runtime_seconds;
    total_runtime = trial.runtime_seconds;
end

function local_print_parameter_comparison(trial)
    truth = trial.true_parameters;
    fprintf('True physical paths:\n');
    fprintf('  Path    Delay    Doppler       AoD       AoA       Gain (real, imag)\n');
    for index = 1:numel(truth.pathDelays)
        fprintf('  %4d  %7.3f  %9.3f  %8.3f  %8.3f  %10.4f %+10.4fj\n', ...
            index, truth.pathDelays(index), truth.pathDopplers(index), ...
            truth.AoDs(index), truth.AoAs(index), ...
            real(truth.pathGains(index)), imag(truth.pathGains(index)));
    end

    estimate = trial.estimated_parameters;
    fprintf('PA-BL strongest estimated DD atoms:\n');
    fprintf('  Atom    Delay    Doppler      Lambda     ||H_q||_F\n');
    for index = 1:numel(estimate.atom_indices)
        fprintf('  %4d  %7.3f  %9.3f  %10.4g  %12.4g\n', ...
            estimate.atom_indices(index), estimate.pathDelays(index), ...
            estimate.pathDopplers(index), estimate.lambda(index), ...
            estimate.coefficient_norm(index));
    end
    fprintf('  AoA/AoD are not estimated by PA-BL; they are shown only for the true paths.\n\n');
end

function value = local_verbose(params)
    value = isfield(params, 'verbose_output') && params.verbose_output;
end

function local_update_progress(action, total)
    persistent count start_time
    if strcmp(action, 'reset')
        count = 0;
        start_time = tic;
        return;
    end
    count = count + 1;
    fprintf('\rProgress: %d/%d (%.1f s)', count, total, toc(start_time));
end
