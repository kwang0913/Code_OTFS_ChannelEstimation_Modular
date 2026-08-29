function results = run_sim_v2(params)
%RUN_SIM_V2 Monte Carlo simulation loop (SNR x trials).
%   results = RUN_SIM_V2(params) runs the SNR x trials simulation and returns
%   a struct with fields: NMSE_before, NMSE_after, BER_after,
%   SNR_db_list, params.

    use_parallel = params.trials > 1;
    % Monte Carlo RNG seeds. Default preserves the historical 1:trials sequence.
    if ~isfield(params, 'rng_seed_start')
        params.rng_seed_start = 1;
    end
    rng_seeds = (params.rng_seed_start:params.rng_seed_start + params.trials - 1).';

    % Initialize results storage (coarse vs final estimates)
    nSNR = length(params.SNR_db_list);
    NMSE_before = NaN(params.trials, nSNR);
    NMSE_after  = NaN(params.trials, nSNR);
    BER_after   = zeros(params.trials, nSNR); % zero sentinel when multi-trial BER is skipped
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
        temp_ber_after   = zeros(params.trials, 1);
        temp_pct_coarse  = NaN(params.trials, K, 5);
        temp_pct_ssr     = NaN(params.trials, K, 5);
        temp_n_detected  = zeros(params.trials, 1);
        temp_coarse_ests = NaN(params.trials, K, 5);

        if use_parallel
            dq = parallel.pool.DataQueue;
            afterEach(dq, @(~) update_progress());
            params.show_figures = false;
            params.verbose_output = false;
            % Capture main session RNG state so all workers start from the same point,
            % ensuring rng(trial_idx) produces identical sequences in parfor and serial.
            mainRngState = rng();
            parfor (trial_idx = 1:params.trials)
                rng(mainRngState);  % reset worker RNG to main session state
                rng(rng_seeds(trial_idx));
                [chanParams, derivedParams] = generate_channel_params(params);
                [XTF, ~, params_t] = generate_TF(params);
                [YTF, ~, ~] = channel_and_noise_rect(XTF, chanParams, params_t, derivedParams, SNR_db);
                coarse_ests = estimate_coarse_params_v2(YTF, XTF, chanParams, params_t, n0, derivedParams);
                final_ests = estimate_ssr_params_constrained(YTF, XTF, coarse_ests, params_t, n0, chanParams);
                temp_nmse_before(trial_idx) = compute_nmse(chanParams, coarse_ests, params_t, derivedParams);
                temp_nmse_after(trial_idx)  = compute_nmse(chanParams, final_ests, params_t, derivedParams);
                % Calculate how good is coarse estimation
                [temp_pct_coarse(trial_idx, :, :), temp_coarse_ests(trial_idx, :, :)] = ...
                    compute_pct_errors(chanParams, coarse_ests);
                temp_pct_ssr(trial_idx, :, :) = compute_pct_errors(chanParams, final_ests);
                temp_n_detected(trial_idx) = coarse_ests.n_detected;
                send(dq, 0);
            end
        else
            for trial_idx = 1:params.trials
                rng(rng_seeds(trial_idx));
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
                temp_ber_after(trial_idx)  = detect_lmmse(rx_OTFS, XTF, params_t, final_ests, derivedParams, databits, guard_mask, n0);

                if params.verbose_output
                    fprintf('  Trial %d: NMSE=%.2f→%.2f dB, BER=%.4e\n', ...
                        trial_idx, ...
                        10*log10(temp_nmse_before(trial_idx)), ...
                        10*log10(temp_nmse_after(trial_idx)), ...
                        temp_ber_after(trial_idx));
                end
            end
        end
        fprintf('\n');

        NMSE_before(:, snr_idx) = temp_nmse_before;
        NMSE_after(:, snr_idx)  = temp_nmse_after;
        BER_after(:, snr_idx)   = temp_ber_after;
        n_detected_aoa(:, snr_idx) = temp_n_detected;
        for ki = 1:K
            for pi = 1:5
                pct_err_coarse(:, snr_idx, ki, pi) = temp_pct_coarse(:, ki, pi);
                pct_err_ssr(:, snr_idx, ki, pi)    = temp_pct_ssr(:, ki, pi);
            end
        end
        for ki = 1:K
            for pi = 1:5
                coarse_ests_all(:, snr_idx, ki, pi) = temp_coarse_ests(:, ki, pi);
            end
        end
        fprintf('BER  : %.4e\n', mean(temp_ber_after));
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

function update_progress()
    persistent count t0
    if isempty(count)
        count = 0;
        t0 = tic;
    end
    count = count + 1;
    fprintf('\rProgress: %d (%.1f s)', count, toc(t0));
end
