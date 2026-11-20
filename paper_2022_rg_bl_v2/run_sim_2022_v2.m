function results = run_sim_2022_v2(params)
%RUN_SIM_2022_V2 Aggregate main-physical Srivastava 2022 trials.
    local_validate(params);
    % Per-trial seeding switches the generator; restore the caller's state.
    caller_rng = rng;
    restore_rng = onCleanup(@() rng(caller_rng)); %#ok<NASGU>
    % compute_ber empty: detect only for a single-trial run.
    if isempty(params.compute_ber), params.compute_ber = params.trials == 1; end
    validateattributes(params.compute_ber, {'logical','numeric'}, {'scalar'}, mfilename, 'compute_ber');
    params.compute_ber = logical(params.compute_ber);
    n_trials = params.trials;
    n_snr = numel(params.SNR_db_list);
    names = {'NMSE_RGBL','NMSE_RGOMP','coefficient_NMSE_RGBL', ...
        'coefficient_NMSE_RGOMP','BER_RGBL','BER_RGOMP', ...
        'support_precision_RGBL','support_recall_RGBL','support_f1_RGBL', ...
        'support_precision_RGOMP','support_recall_RGOMP','support_f1_RGOMP', ...
        'rgbl_runtime_seconds','rgomp_runtime_seconds', ...
        'detection_runtime_seconds','runtime_seconds'};
    for index = 1:numel(names)
        metrics.(names{index}) = zeros(n_trials,n_snr);
    end
    metrics.rgbl_iterations = zeros(n_trials,n_snr);
    metrics.rgomp_iterations = zeros(n_trials,n_snr);
    metrics.rgbl_converged = false(n_trials,n_snr);
    metrics.rgomp_converged = false(n_trials,n_snr);
    trials = cell(n_trials,n_snr);
    parallel_available = local_parallel_available();
    use_parallel = params.use_parallel && n_trials > 1 && parallel_available;
    trial_params = params;
    trial_params.verbose_output = false;
    for snr_index = 1:n_snr
        snr_db = params.SNR_db_list(snr_index);
        fprintf('--- Srivastava 2022 v2 SNR = %g dB ---\n',snr_db);
        temporary = cell(n_trials,1);
        local_update_progress('reset',n_trials);
        if use_parallel
            dq = parallel.pool.DataQueue;
            afterEach(dq,@(~) local_update_progress('increment',n_trials));
            parfor trial_index = 1:n_trials
                temporary{trial_index} = simulate_trial_2022_v2( ...
                    trial_params,snr_db,trial_index);
                send(dq,0);
            end
        else
            for trial_index = 1:n_trials
                temporary{trial_index} = simulate_trial_2022_v2( ...
                    trial_params,snr_db,trial_index);
                local_update_progress('increment',n_trials);
                if local_verbose(params)
                    trial = temporary{trial_index};
                    fprintf(['\n  Trial %d: RG-BL NMSE=%.2f dB, ' ...
                        'RG-OMP NMSE=%.2f dB, BER=%.4e / %.4e\n'], ...
                        trial_index,10*log10(trial.NMSE_RGBL), ...
                        10*log10(trial.NMSE_RGOMP), ...
                        trial.BER_RGBL,trial.BER_RGOMP);
                end
            end
        end
        fprintf('\n');
        for trial_index = 1:n_trials
            trial = temporary{trial_index};
            direct = {'NMSE_RGBL','NMSE_RGOMP','coefficient_NMSE_RGBL', ...
                'coefficient_NMSE_RGOMP','BER_RGBL','BER_RGOMP', ...
                'rgbl_runtime_seconds','rgomp_runtime_seconds', ...
                'detection_runtime_seconds','runtime_seconds'};
            for index = 1:numel(direct)
                metrics.(direct{index})(trial_index,snr_index) = ...
                    trial.(direct{index});
            end
            metrics.support_precision_RGBL(trial_index,snr_index) = ...
                trial.support_RGBL.precision;
            metrics.support_recall_RGBL(trial_index,snr_index) = ...
                trial.support_RGBL.recall;
            metrics.support_f1_RGBL(trial_index,snr_index) = ...
                trial.support_RGBL.f1;
            metrics.support_precision_RGOMP(trial_index,snr_index) = ...
                trial.support_RGOMP.precision;
            metrics.support_recall_RGOMP(trial_index,snr_index) = ...
                trial.support_RGOMP.recall;
            metrics.support_f1_RGOMP(trial_index,snr_index) = ...
                trial.support_RGOMP.f1;
            metrics.rgbl_iterations(trial_index,snr_index) = ...
                trial.rgbl_state.iterations;
            metrics.rgomp_iterations(trial_index,snr_index) = ...
                trial.rgomp_state.iterations;
            metrics.rgbl_converged(trial_index,snr_index) = ...
                trial.rgbl_state.converged;
            metrics.rgomp_converged(trial_index,snr_index) = ...
                trial.rgomp_state.converged;
        end
        trials(:,snr_index) = temporary;
        if n_trials == 1
            local_print_parameter_comparison(temporary{1});
        end
        fprintf('Main NMSE (RG-BL / RG-OMP): %.2f / %.2f dB\n', ...
            10*log10(mean(metrics.NMSE_RGBL(:,snr_index))), ...
            10*log10(mean(metrics.NMSE_RGOMP(:,snr_index))));
        if params.compute_ber
            fprintf('BER (RG-BL / RG-OMP): %.4e / %.4e\n\n', ...
                mean(metrics.BER_RGBL(:,snr_index)), ...
                mean(metrics.BER_RGOMP(:,snr_index)));
        else
            fprintf('BER (RG-BL / RG-OMP) skipped (compute_ber = false).\n\n');
        end
    end
    results = metrics;
    results.trials = trials;
    results.SNR_db_list = params.SNR_db_list;
    results.params = params;
end

function local_validate(params)
    if ~isstruct(params) || ~isscalar(params) || params.trials < 1 || ...
            params.trials ~= round(params.trials) || isempty(params.SNR_db_list) || ...
            any(~isfinite(params.SNR_db_list(:)))
        error('run_sim_2022_v2:InvalidParams','Invalid trials or SNR_db_list.');
    end
end

function available = local_parallel_available()
    available = false;
    try
        available = license('test','Parallel_Computing_Toolbox') || ...
            license('test','Distrib_Computing_Toolbox');
    catch
    end
end

function local_update_progress(action,total)
    persistent count start_time
    if strcmp(action,'reset')
        count = 0;
        start_time = tic;
        return;
    end
    count = count+1;
    fprintf('\rProgress: %d/%d (%.1f s)',count,total,toc(start_time));
end

function local_print_parameter_comparison(trial)
    truth = trial.true_parameters;
    fprintf('True physical paths:\n');
    fprintf('  Path    Delay    Doppler       AoD       AoA       Gain (real, imag)\n');
    for index = 1:numel(truth.pathDelays)
        fprintf('  %4d  %7.3f  %9.3f  %8.3f  %8.3f  %10.4f %+10.4fj\n', ...
            index,truth.pathDelays(index),truth.pathDopplers(index), ...
            truth.AoDs(index),truth.AoAs(index), ...
            real(truth.pathGains(index)),imag(truth.pathGains(index)));
    end
    local_print_estimate('RG-BL',trial.estimated_parameters.RGBL);
    local_print_estimate('RG-OMP',trial.estimated_parameters.RGOMP);
    fprintf('  AoA/AoD are not estimated; they are shown only for the true paths.\n\n');
end

function local_print_estimate(name,estimate)
    fprintf('%s strongest estimated DD atoms:\n',name);
    fprintf('  Atom    Delay    Doppler       Score     ||H_q||_F\n');
    for index = 1:numel(estimate.atom_indices)
        fprintf('  %4d  %7.3f  %9.3f  %10.4g  %12.4g\n', ...
            estimate.atom_indices(index),estimate.pathDelays(index), ...
            estimate.pathDopplers(index),estimate.score(index), ...
            estimate.coefficient_norm(index));
    end
end

function value = local_verbose(params)
    value = isfield(params,'verbose_output') && params.verbose_output;
end
