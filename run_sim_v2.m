function results = run_sim_v2(params)
%RUN_SIM_V2 Monte Carlo simulation loop (SNR x trials).
%   results = RUN_SIM_V2(params) runs the SNR x trials simulation and returns
%   a struct with fields: NMSE_before, NMSE_after, BER_after,
%   SNR_db_list, params.
%
%   params.use_parallel runs trials with parfor when trials > 1. Trial i is
%   seeded with rng(rng_seed_start + i - 1, 'threefry') in both modes.
%   params.compute_ber (see should_compute_ber) controls data detection;
%   skipped trials report NaN BER.

    % Per-trial seeding switches the generator; restore the caller's state.
    caller_rng = rng;
    restore_rng = onCleanup(@() rng(caller_rng)); %#ok<NASGU>
    params.compute_ber = should_compute_ber(params);
    compute_ber = params.compute_ber;
    use_parallel = params.use_parallel && params.trials > 1;
    if use_parallel && ~local_parallel_available()
        warning('run_sim_v2:ParallelUnavailable', ...
            'Parallel execution was requested but the Parallel Computing Toolbox is unavailable; using serial execution.');
        use_parallel = false;
    end
    rng_seeds = (params.rng_seed_start:params.rng_seed_start + params.trials - 1).';

    % Initialize results storage (coarse vs final estimates)
    nSNR = length(params.SNR_db_list);
    NMSE_before = NaN(params.trials, nSNR);
    NMSE_after  = NaN(params.trials, nSNR);
    BER_after   = NaN(params.trials, nSNR); % NaN when detection is skipped
    K = params.targets_num;
    % Pct errors: (trials, nSNR, K, 5) where dim4 = [AoA, AoD, delay, doppler, gain]
    pct_err_coarse = NaN(params.trials, nSNR, K, 5);
    pct_err_ssr    = NaN(params.trials, nSNR, K, 5);
    n_detected_aoa = zeros(params.trials, nSNR);
    % Per-target coarse estimates: (trials, nSNR, K, 5) = [AoA, AoD, delay, doppler, gain]
    coarse_ests_all = NaN(params.trials, nSNR, K, 5);

    for snr_idx = 1:length(params.SNR_db_list)
        SNR_db = params.SNR_db_list(snr_idx);
        n0 = 1 / (10^(SNR_db / 10));
        fprintf('--- SNR = %d dB ---\n', SNR_db);

        temp_nmse_before = NaN(params.trials, 1);
        temp_nmse_after  = NaN(params.trials, 1);
        temp_ber_after   = NaN(params.trials, 1);
        temp_pct_coarse  = NaN(params.trials, K, 5);
        temp_pct_ssr     = NaN(params.trials, K, 5);
        temp_n_detected  = zeros(params.trials, 1);
        temp_coarse_ests = NaN(params.trials, K, 5);

        if use_parallel
            update_progress('reset', params.trials);
            dq = parallel.pool.DataQueue;
            afterEach(dq, @(~) update_progress('increment', params.trials));
            params.show_figures = false;
            params.verbose_output = false;
            parfor (trial_idx = 1:params.trials)
                rng(rng_seeds(trial_idx), 'threefry');
                [chanParams, derivedParams] = generate_channel_params(params);
                [XTF, ~, params_t, databits, guard_mask] = generate_TF(params);
                [YTF, ~, rx_OTFS] = channel_and_noise_rect(XTF, chanParams, params_t, derivedParams, SNR_db);
                coarse_ests = estimate_coarse_params_v2(YTF, XTF, chanParams, params_t, n0, derivedParams);
                final_ests = estimate_ssr_params_constrained(YTF, XTF, coarse_ests, params_t, n0, chanParams);
                temp_nmse_before(trial_idx) = compute_nmse(chanParams, coarse_ests, params_t, derivedParams);
                temp_nmse_after(trial_idx)  = compute_nmse(chanParams, final_ests, params_t, derivedParams);
                % Calculate how good is coarse estimation
                [temp_pct_coarse(trial_idx, :, :), temp_coarse_ests(trial_idx, :, :)] = ...
                    compute_pct_errors(chanParams, coarse_ests);
                temp_pct_ssr(trial_idx, :, :) = compute_pct_errors(chanParams, final_ests);
                temp_n_detected(trial_idx) = coarse_ests.n_detected;
                if compute_ber
                    temp_ber_after(trial_idx) = detect_lmmse(rx_OTFS, XTF, params_t, final_ests, derivedParams, databits, guard_mask, n0);
                end
                send(dq, 0);
            end
        else
            for trial_idx = 1:params.trials
                rng(rng_seeds(trial_idx), 'threefry');
                [chanParams, derivedParams] = generate_channel_params(params);
                [XTF, ~, params_t, databits, guard_mask] = generate_TF(params);
                [YTF, ~, rx_OTFS] = channel_and_noise_rect(XTF, chanParams, params_t, derivedParams, SNR_db);
                coarse_ests = estimate_coarse_params_v2(YTF, XTF, chanParams, params_t, n0, derivedParams);
                final_ests = estimate_ssr_params_constrained(YTF, XTF, coarse_ests, params_t, n0, chanParams);
                temp_nmse_before(trial_idx) = compute_nmse(chanParams, coarse_ests, params_t, derivedParams);
                temp_nmse_after(trial_idx)  = compute_nmse(chanParams, final_ests, params_t, derivedParams);
                [temp_pct_coarse(trial_idx, :, :), temp_coarse_ests(trial_idx, :, :)] = ...
                    compute_pct_errors(chanParams, coarse_ests);
                temp_pct_ssr(trial_idx, :, :) = compute_pct_errors(chanParams, final_ests);
                temp_n_detected(trial_idx) = coarse_ests.n_detected;
                if compute_ber
                    temp_ber_after(trial_idx) = detect_lmmse(rx_OTFS, XTF, params_t, final_ests, derivedParams, databits, guard_mask, n0);
                end

                if params.verbose_output
                    fprintf('  Trial %d: NMSE=%.2f→%.2f dB', trial_idx, ...
                        10*log10(temp_nmse_before(trial_idx)), ...
                        10*log10(temp_nmse_after(trial_idx)));
                    if compute_ber
                        fprintf(', BER=%.4e', temp_ber_after(trial_idx));
                    end
                    fprintf('\n');
                end
            end
        end
        fprintf('\n');

        NMSE_before(:, snr_idx) = temp_nmse_before;
        NMSE_after(:, snr_idx)  = temp_nmse_after;
        BER_after(:, snr_idx)   = temp_ber_after;
        n_detected_aoa(:, snr_idx) = temp_n_detected;
        pct_err_coarse(:, snr_idx, :, :)  = reshape(temp_pct_coarse, params.trials, 1, K, 5);
        pct_err_ssr(:, snr_idx, :, :)     = reshape(temp_pct_ssr, params.trials, 1, K, 5);
        coarse_ests_all(:, snr_idx, :, :) = reshape(temp_coarse_ests, params.trials, 1, K, 5);
        if compute_ber
            fprintf('BER  : %.4e\n', mean(temp_ber_after));
        else
            fprintf('BER  : skipped (compute_ber = false)\n');
        end
        fprintf('NMSE : %.2f -> %.2f dB\n', ...
            10*log10(mean(temp_nmse_before, 'omitnan')), ...
            10*log10(mean(temp_nmse_after, 'omitnan')));
        fprintf('\n');
    end

    % Package results
    results.NMSE_before    = NMSE_before;
    results.NMSE_after     = NMSE_after;
    results.BER_after      = BER_after;
    results.pct_err_coarse = pct_err_coarse;  % (trials, nSNR, K, 5): [AoA, AoD, delay, doppler, gain]
    results.pct_err_ssr    = pct_err_ssr;
    results.n_detected_aoa = n_detected_aoa;  % (trials, nSNR)
    results.coarse_ests_all = coarse_ests_all;  % (trials, nSNR, K, 5): [AoA, AoD, delay, doppler, gain]
    results.SNR_db_list    = params.SNR_db_list;
    results.params = params;
    results.rng_seeds = rng_seeds;
end

function available = local_parallel_available()
    available = false;
    try
        available = license('test', 'Parallel_Computing_Toolbox') || ...
            license('test', 'Distrib_Computing_Toolbox');
    catch
    end
end

function update_progress(action, total)
    persistent count t0
    if strcmp(action, 'reset')
        count = 0;
        t0 = tic;
        return;
    end
    count = count + 1;
    fprintf('\rProgress: %d/%d (%.1f s)', count, total, toc(t0));
end
