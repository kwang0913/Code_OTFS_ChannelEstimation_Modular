% run_comparison_nr1_v2.m
% Nt=16 TX, Nr=1 RX. L-shaped pilots + DFT CLEAN coarse + block-OMP SSR.
% Supports 3GPP clustered channels (Ns subpaths per cluster) or point targets (Ns=1).
% Oracle noise as in B (actual receive-noise variance). LMMSE detection (BER) runs only
% for a single-trial run (trials = 1).

close all;
clearvars;

% --- Configuration ---
ov = struct();
ov.Nt = 16;
ov.Nr = 1;
ov.N  = 140;
ov.df = 15e3;
ov.trials = 200;
ov.SNR_db_list = 0:2:20;
ov.pilot_power = 8;
ov.selected_targets = [1,2,3,4];
ov.target_aods = [-31, 46, -46, 20];  % degrees, one per target
ov.n_targets = [];  % fixed pool targets selected_targets (not a random draw per trial)
ov.Ns_subpaths = 20;
ov.angular_spread_deg = 2;

params = config_params_std(ov);

% L-shaped pilot layout: each TX gets a delay arm + Doppler arm sharing a corner.
% Both arms per TX resolve delay AND Doppler; TX orthogonality via distinct corners.
b = 121;
Nt = params.Nt;
b_arm = floor(b / 2);   % pilots per arm (60 for b=121)

% Strides: spread each arm across the full dimension
m_stride = max(1, floor((params.M - 1) / max(b_arm - 1, 1)));
n_stride = max(1, floor((params.N - 1) / max(b_arm - 1, 1)));

% Corner positions: spread Nt corners across (M, N) for diversity
m_corners = round(linspace(1, params.M - (b_arm-1)*m_stride, Nt));
n_corners = round(linspace(1, params.N - (b_arm-1)*n_stride, Nt));

pb = [];
blk = 0;
for tx = 1:Nt
    mc = m_corners(tx);
    nc = n_corners(tx);

    % Block A: delay arm — b_arm m-positions at fixed n = nc
    blk = blk + 1;
    pb(blk).dim    = 'delay';
    pb(blk).m      = mc + (0:b_arm-1) * m_stride;
    pb(blk).n      = nc;
    pb(blk).tx     = tx;
    pb(blk).ref_m  = mc;
    pb(blk).ref_n  = nc;
    pb(blk).ref_tx = tx;

    % Block B: Doppler arm — b_arm n-positions at fixed m = mc
    blk = blk + 1;
    pb(blk).dim    = 'doppler';
    pb(blk).m      = mc;
    pb(blk).n      = nc + (0:b_arm-1) * n_stride;
    pb(blk).tx     = tx;
    pb(blk).ref_m  = mc;
    pb(blk).ref_n  = nc;
    pb(blk).ref_tx = tx;
end

params.pilot_blocks_def = pb;  % config_pilots derives offsets and pilot counts
params.n_meas = compute_n_meas(config_pilots(params), params.Nr);

% --- Run Simulation ---
results = run_sim_nr1(params);

% --- Save Results ---
if params.trials > 1
    out_dir = 'results';
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    fname = sprintf('comparison_nr1_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
    save(fullfile(out_dir, fname), 'results');
    fprintf('Saved: %s\n', fname);
end

% --- Plot NMSE ---
figure;
plot(results.SNR_db_list, 10*log10(mean(results.nmse_coarse, 1)), '-o', 'LineWidth', 2); hold on;
plot(results.SNR_db_list, 10*log10(mean(results.nmse_ssr, 1)),  '-s', 'LineWidth', 2);
grid on;
xlabel('SNR [dB]'); ylabel('NMSE [dB]');
legend('OMP', 'OMP + Refined');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Nt=%d, Nr=1 Comparison', params.Nt));

% --- Plot BER (LMMSE is evaluated only for a single-trial run) ---
if params.trials == 1
    figure;
    ber_plot = results.BER(1,:);
    zero_ber = ber_plot == 0;
    ber_plot(zero_ber) = NaN; % A logarithmic axis cannot represent zero.
    semilogy(results.SNR_db_list, ber_plot, '-s', 'LineWidth', 2);
    grid on;
    xlabel('SNR [dB]'); ylabel('BER');
    legend('OMP + Refined + LMMSE', 'Location', 'best');
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Nt=%d, Nr=1 BER', params.Nt));
    finite_snr = results.SNR_db_list(isfinite(results.SNR_db_list));
    if ~isempty(finite_snr)
        limits = [min(finite_snr), max(finite_snr)];
        if limits(1) == limits(2), limits = limits + [-1, 1]; end
        xlim(limits);
    end
    if any(zero_ber)
        if all(zero_ber), ylim([1e-6, 1]); end
        text(0.03, 0.07, sprintf('BER = 0 at SNR [dB]: %s\n(zero points omitted on log axis)', ...
            strtrim(sprintf('%g ', results.SNR_db_list(zero_ber)))), ...
            'Units', 'normalized', 'VerticalAlignment', 'bottom');
    end
end

% --- Display Per-Target Errors ---
param_labels = {'Angle', 'Delay', 'Doppler'};
K = size(results.err_coarse, 3);
nSNR = length(results.SNR_db_list);

fprintf('\n--- Per-Target Abs Errors (mean +/- std across trials) ---\n');
fprintf('%-16s', 'SNR(dB)');
fprintf('%16d ', results.SNR_db_list);
fprintf('\n');
fprintf('  OMP (init):\n');
for k = 1:K
    fprintf('    Target %d:\n', k);
    for p = [2, 3]
        fprintf('    %-12s', param_labels{p});
        for s = 1:nSNR
            m = mean(results.err_coarse(:, s, k, p), 1, 'omitnan');
            sd = std(results.err_coarse(:, s, k, p), 0, 1, 'omitnan');
            fprintf('%7.2f+/-%-7.2f', m, sd);
        end
        fprintf('\n');
    end
end
fprintf('  OMP (refined):\n');
for k = 1:K
    fprintf('    Target %d:\n', k);
    for p = [2, 3]
        fprintf('    %-12s', param_labels{p});
        for s = 1:nSNR
            m = mean(results.err_ssr(:, s, k, p), 1, 'omitnan');
            sd = std(results.err_ssr(:, s, k, p), 0, 1, 'omitnan');
            fprintf('%7.2f+/-%-7.2f', m, sd);
        end
        fprintf('\n');
    end
end
fprintf('--------------------------------------\n\n');


%% ========================================================================
%  Simulation loop
%  ========================================================================

function results = run_sim_nr1(params)
%RUN_SIM_NR1 SNR x trials loop. Trial i uses rng(rng_seed_start + i - 1, 'threefry')
%   at every SNR; the caller's RNG state is restored on exit. LMMSE detection
%   (BER) runs only when trials == 1; otherwise BER stays NaN.
    caller_rng = rng;
    restore_rng = onCleanup(@() rng(caller_rng)); %#ok<NASGU>
    nSNR = length(params.SNR_db_list);
    K = params.targets_num;
    use_parallel = params.trials > 1 && params.n_workers > 0;
    compute_ber = params.trials == 1;
    rng_seeds = (params.rng_seed_start:params.rng_seed_start + params.trials - 1).';
    derivedParams = compute_derived_params(params);

    nmse_coarse = NaN(params.trials, nSNR);
    nmse_ssr    = NaN(params.trials, nSNR);
    BER         = NaN(params.trials, nSNR); % NaN means detection was skipped.
    err_coarse  = NaN(params.trials, nSNR, K, 3);
    err_ssr     = NaN(params.trials, nSNR, K, 3);

    for snr_idx = 1:nSNR
        SNR_db = params.SNR_db_list(snr_idx);
        fprintf('--- SNR = %d dB ---\n', SNR_db);

        tr = cell(params.trials, 1);
        if use_parallel
            update_progress(params.trials);
            dq = parallel.pool.DataQueue;
            afterEach(dq, @(~) update_progress());
            params_par = params;
            params_par.show_figures = false;
            params_par.verbose_output = false;
            parfor (trial_idx = 1:params.trials, params.n_workers)
                tr{trial_idx} = run_trial_nr1(params_par, derivedParams, SNR_db, rng_seeds(trial_idx), false);
                send(dq, 0);
            end
        else
            for trial_idx = 1:params.trials
                tr{trial_idx} = run_trial_nr1(params, derivedParams, SNR_db, rng_seeds(trial_idx), compute_ber);
            end
        end
        fprintf('\n');

        for trial_idx = 1:params.trials
            r = tr{trial_idx};
            nmse_coarse(trial_idx, snr_idx)      = r.nmse_coarse;
            nmse_ssr(trial_idx, snr_idx)         = r.nmse_ssr;
            BER(trial_idx, snr_idx)              = r.ber;
            err_coarse(trial_idx, snr_idx, :, :) = r.err_coarse;
            err_ssr(trial_idx, snr_idx, :, :)    = r.err_ssr;
        end

        fprintf('NMSE : %.2f -> %.2f dB\n', ...
            10*log10(mean(nmse_coarse(:, snr_idx))), 10*log10(mean(nmse_ssr(:, snr_idx))));
        if compute_ber
            fprintf('BER  : %.4e\n', BER(1, snr_idx));
        end
        fprintf('\n');
    end

    results.nmse_coarse = nmse_coarse;
    results.nmse_ssr    = nmse_ssr;
    results.BER         = BER;
    results.err_coarse  = err_coarse;  % (trials, nSNR, K, 3): [angle, delay, doppler]
    results.err_ssr     = err_ssr;
    results.rng_seeds   = rng_seeds;
    results.SNR_db_list = params.SNR_db_list;
    results.params      = params;
end

function r = run_trial_nr1(params, derivedParams, SNR_db, seed, compute_ber)
%RUN_TRIAL_NR1 One seeded trial: clustered channel -> DFT CLEAN -> block-OMP SSR -> NMSE and errors.
%   COMPUTE_BER adds LMMSE detection with the SSR estimate.
%   As in B, this Nr = 1 comparison is oracle-noise: estimators and detector get
%   the actual receive-noise variance n0 (run_sim uses the estimate from rx_OTFS).
    rng(seed, 'threefry');
    chanParams = generate_channel_params_nr1(params);
    [XTF, ~, params_t, databits, guard_mask] = generate_TF(params);
    % n0 is the actual receive-noise variance, before any pilot normalization.
    [YTF, ~, rx_OTFS, n0] = channel_and_noise_rect(XTF, chanParams, params_t, derivedParams, SNR_db);

    [coarse_ests, profiles] = estimate_coarse_dft_nr1(YTF, XTF, params_t, n0);

    if params.show_figures
        % True cluster means for plotting
        true_del_plot = chanParams.mean_delays;
        true_dop_plot = chanParams.mean_dopplers;

        % --- Delay CLEAN iterations (2x2 subplot) ---
        n_di = profiles.n_delay_iters;
        n_sub = min(n_di, 4);
        xl_d = [0 max(true_del_plot) + 20];
        figure;
        for si = 1:n_sub
            subplot(2, 2, si); hold on;
            plot(profiles.delay_axis, profiles.delay_clean_hist(si,:), '-', 'LineWidth', 1.5);
            for jj = 1:length(true_del_plot)
                xline(true_del_plot(jj), 'r--', 'LineWidth', 1.5);
            end
            if si <= length(profiles.cand_delays)
                xline(profiles.cand_delays(si), 'g:', 'LineWidth', 2);
            end
            hold off; xlim(xl_d); grid on;
            title(sprintf('Delay CLEAN iter %d', si));
            xlabel('Delay [tap]'); ylabel('Norm. Mag.');
            set(gca, 'FontSize', 11, 'FontWeight', 'bold');
        end
        sgtitle('Delay CLEAN Successive Cancellation', 'FontSize', 13, 'FontWeight', 'bold');

        % --- Doppler CLEAN iterations (2x2 subplot) ---
        n_vi = profiles.n_doppler_iters;
        n_sub = min(n_vi, 4);
        figure;
        for si = 1:n_sub
            subplot(2, 2, si); hold on;
            plot(profiles.doppler_axis, profiles.doppler_clean_hist(si,:), '-', 'LineWidth', 1.5);
            for jj = 1:length(true_dop_plot)
                xline(true_dop_plot(jj), 'r--', 'LineWidth', 1.5);
            end
            if si <= length(profiles.cand_dopplers)
                xline(profiles.cand_dopplers(si), 'g:', 'LineWidth', 2);
            end
            hold off; grid on;
            title(sprintf('Doppler CLEAN iter %d', si));
            xlabel('Doppler [tap]'); ylabel('Norm. Mag.');
            set(gca, 'FontSize', 11, 'FontWeight', 'bold');
        end
        sgtitle('Doppler CLEAN Successive Cancellation', 'FontSize', 13, 'FontWeight', 'bold');
    end

    [final_ests, init_ests] = estimate_ssr_nr1(YTF, XTF, params_t, derivedParams, n0, coarse_ests);

    init_ests.beta_per_tx = per_antenna_ls(YTF, XTF, init_ests, params_t, derivedParams, n0);
    r.nmse_coarse = compute_nmse_nr1(chanParams, init_ests, params_t, derivedParams);
    r.nmse_ssr    = compute_nmse_nr1(chanParams, final_ests, params_t, derivedParams);
    truth = cluster_means_nr1(chanParams);
    r.err_coarse  = compute_errors(truth, init_ests, params_t);
    r.err_ssr     = compute_errors(truth, final_ests, params_t);
    r.err_coarse  = r.err_coarse(:, 1:3);
    r.err_ssr     = r.err_ssr(:, 1:3);

    r.ber = NaN;
    if compute_ber
        r.ber = detect_lmmse(rx_OTFS, XTF, params_t, per_tx_channel_nr1(final_ests), derivedParams, ...
            databits, guard_mask, [], n0);
    end

    if params.verbose_output
        fprintf('  NMSE=%.2f -> %.2f dB\n', 10*log10(r.nmse_coarse), 10*log10(r.nmse_ssr));
    end
end

function update_progress(total)
%UPDATE_PROGRESS update_progress(total) resets the counter; update_progress() increments.
    persistent count n_total t0
    if nargin
        count = 0;  n_total = total;  t0 = tic;
        return
    end
    count = count + 1;
    fprintf('\rProgress: %d/%d (%.1f s)', count, n_total, toc(t0));
end


%% ========================================================================
%  Coarse estimation: DFT CLEAN on L-shaped pilot arms
%  ========================================================================

function [coarse_ests, profiles] = estimate_coarse_dft_nr1(YTF, XTF, params, n0)
%ESTIMATE_COARSE_DFT_NR1 DFT-based coarse delay+Doppler estimation.
%   Delay CLEAN on all delay arms, Doppler CLEAN on all Doppler arms.
%   Returns unpaired candidates for SSR block-OMP to pair jointly.

    K  = params.targets_num;
    M  = params.M;
    N  = params.N;
    pb = params.pilot_blocks;
    Y_single = squeeze(YTF(:,:,1));

    P = 2048;

    % Identify delay arms (many m, single n) and Doppler arms (single m, many n)
    delay_blk_idx   = [];
    doppler_blk_idx = [];
    for bi = 1:length(pb)
        if length(pb(bi).m) > 1 && length(pb(bi).n) == 1
            delay_blk_idx(end+1) = bi; %#ok<AGROW>
        elseif length(pb(bi).m) == 1 && length(pb(bi).n) > 1
            doppler_blk_idx(end+1) = bi; %#ok<AGROW>
        end
    end
    n_delay_arms   = length(delay_blk_idx);
    n_doppler_arms = length(doppler_blk_idx);

    %% 1. Delay CLEAN: IFFT each arm, average spectra, successive cancellation
    b_d = length(pb(delay_blk_idx(1)).m);
    stride_d = 1;
    if b_d > 1, stride_d = pb(delay_blk_idx(1)).m(2) - pb(delay_blk_idx(1)).m(1); end

    delay_axis  = (0:P-1) * M / (P * stride_d);  % physical delay axis (taps)
    max_unamb_d = M / stride_d;                   % unambiguous range set by pilot spacing
    valid_d     = delay_axis < max_unamb_d;

    delay_ratios = zeros(b_d, n_delay_arms);
    for i = 1:n_delay_arms
        blk = pb(delay_blk_idx(i));
        delay_ratios(:,i) = Y_single(blk.m, blk.n(1)) ./ XTF(blk.m, blk.n(1), blk.tx);
    end
    delay_ratios_clean = delay_ratios;

    cand_delays      = [];
    orig_peak_d      = 0;
    avg_delay_init   = [];
    delay_clean_hist = zeros(K, P);
    n_delay_iters    = 0;

    for iter = 1:K
        specs = zeros(n_delay_arms, P);
        for i = 1:n_delay_arms
            specs(i,:) = abs(ifft(delay_ratios_clean(:,i), P));
        end
        avg_spec = mean(specs, 1);
        avg_spec(~valid_d) = 0;

        if iter == 1
            orig_peak_d = max(avg_spec);
            avg_delay_init = avg_spec / orig_peak_d;
        end
        delay_clean_hist(iter,:) = avg_spec / orig_peak_d;

        [pks, locs, ~, proms] = findpeaks(avg_spec, 'SortStr', 'descend');
        if isempty(pks), break; end
        if proms(1) / orig_peak_d < 0.1 && iter > 1, break; end

        spec_norm = avg_spec / orig_peak_d;
        spec_mean = mean(spec_norm(valid_d));
        spec_med  = median(spec_norm(valid_d));
        if spec_mean < spec_med && ...
                (pks(1)/orig_peak_d - spec_med) < (spec_med - spec_mean) && iter > 1
            break;
        end

        n_delay_iters = iter;
        peak_delay = delay_axis(locs(1));
        cand_delays(end+1) = peak_delay; %#ok<AGROW>

        % LS-subtract detected delay component from all arms (CLEAN cancellation)
        a_d = exp(-1j * 2 * pi * (0:b_d-1).' * stride_d * peak_delay / M);
        for i = 1:n_delay_arms
            c_i = (a_d' * a_d) \ (a_d' * delay_ratios_clean(:,i));
            delay_ratios_clean(:,i) = delay_ratios_clean(:,i) - a_d * c_i;
        end
    end

    if isempty(cand_delays)
        [~, fb] = max(avg_delay_init);
        cand_delays = delay_axis(fb);
    end
    cand_delays = round(cand_delays);

    %% 2. Doppler CLEAN: FFT each arm, average spectra, successive cancellation
    b_v = length(pb(doppler_blk_idx(1)).n);
    stride_v = 1;
    if b_v > 1, stride_v = pb(doppler_blk_idx(1)).n(2) - pb(doppler_blk_idx(1)).n(1); end

    doppler_axis = ((0:P-1) - P/2) * N / (P * stride_v);

    doppler_ratios = zeros(b_v, n_doppler_arms);
    for j = 1:n_doppler_arms
        blk = pb(doppler_blk_idx(j));
        for p = 1:b_v
            doppler_ratios(p,j) = Y_single(blk.m(1), blk.n(p)) / XTF(blk.m(1), blk.n(p), blk.tx);
        end
    end
    doppler_ratios_clean = doppler_ratios;

    cand_dopplers      = [];
    orig_peak_v        = 0;
    avg_doppler_init   = [];
    doppler_clean_hist = zeros(K, P);
    n_doppler_iters    = 0;

    for iter = 1:K
        specs = zeros(n_doppler_arms, P);
        for j = 1:n_doppler_arms
            specs(j,:) = abs(fftshift(fft(doppler_ratios_clean(:,j), P)));
        end
        avg_spec = mean(specs, 1);

        if iter == 1
            orig_peak_v = max(avg_spec);
            avg_doppler_init = avg_spec / orig_peak_v;
        end
        doppler_clean_hist(iter,:) = avg_spec / orig_peak_v;

        [pks, locs, ~, proms] = findpeaks(avg_spec, 'SortStr', 'descend');
        if isempty(pks), break; end
        if proms(1) / orig_peak_v < 0.1 && iter > 1, break; end

        spec_norm = avg_spec / orig_peak_v;
        spec_mean = mean(spec_norm);
        spec_med  = median(spec_norm);
        if spec_mean < spec_med && ...
                (pks(1)/orig_peak_v - spec_med) < (spec_med - spec_mean) && iter > 1
            break;
        end

        n_doppler_iters = iter;
        peak_doppler = doppler_axis(locs(1));
        cand_dopplers(end+1) = peak_doppler; %#ok<AGROW>

        a_v = exp(1j * 2 * pi * (0:b_v-1).' * stride_v * peak_doppler / N);
        for j = 1:n_doppler_arms
            c_j = (a_v' * a_v) \ (a_v' * doppler_ratios_clean(:,j));
            doppler_ratios_clean(:,j) = doppler_ratios_clean(:,j) - a_v * c_j;
        end
    end

    if isempty(cand_dopplers)
        [~, fb] = max(avg_doppler_init);
        cand_dopplers = doppler_axis(fb);
    end

    %% 3. Return unpaired candidates
    coarse_ests = struct('cand_delays', cand_delays, ...
                         'cand_dopplers', cand_dopplers);

    %% 4. Profiles for plotting
    profiles.delay_axis         = delay_axis;
    profiles.delay_profile      = avg_delay_init;
    profiles.delay_clean_hist   = delay_clean_hist(1:max(n_delay_iters,1), :);
    profiles.n_delay_iters      = n_delay_iters;
    profiles.doppler_axis       = doppler_axis;
    profiles.doppler_profile    = avg_doppler_init;
    profiles.doppler_clean_hist = doppler_clean_hist(1:max(n_doppler_iters,1), :);
    profiles.n_doppler_iters    = n_doppler_iters;
    profiles.cand_delays        = cand_delays;
    profiles.cand_dopplers      = cand_dopplers;
    profiles.valid_d            = valid_d;

    if params.verbose_output
        fprintf('--- DFT Coarse (unpaired) ---\n');
        fprintf('Cand. delays:   %s\n', sprintf('%5.1f ', cand_delays));
        fprintf('Cand. dopplers: %s\n', sprintf('%5.2f ', cand_dopplers));
        fprintf('---------------------------------\n');
    end
end


%% ========================================================================
%  SSR: block-OMP on CLEAN-seeded dictionary + binary refinement
%  ========================================================================

function [final_ests, init_ests] = estimate_ssr_nr1(YTF, XTF, params, derivedParams, n0, coarse_ests)
%ESTIMATE_SSR_NR1 Block-OMP with per-antenna coefficients + binary refinement.
%   2D dictionary (delay × Doppler, AoD=0). Per-antenna LS handles AoD.
%   Blocks grouped by TX antenna for joint delay-Doppler resolution.

    K  = params.targets_num;
    N  = params.N;
    Nt = params.Nt;
    pilot_blocks = params.pilot_blocks;

    y = build_measurement_vector(YTF, XTF, pilot_blocks, params.Nr, n0);

    % Group measurement rows (one per unique pilot site, Nr = 1) into per-antenna super-blocks
    sites = pilot_sites(pilot_blocks);
    blk_rows = cell(Nt, 1);
    for tx = 1:Nt
        blk_rows{tx} = find(sites.tx == tx).';
    end
    n_blocks = Nt;

    if nargin >= 6 && ~isempty(coarse_ests) && isfield(coarse_ests, 'cand_delays')
        %% === CLEAN candidates: build dense 2D dictionary ===
        dop_step  = 0.1;
        dop_range = 2;

        all_dops = [];
        for ci = 1:length(coarse_ests.cand_dopplers)
            cv = coarse_ests.cand_dopplers(ci);
            all_dops = [all_dops, (cv - dop_range):dop_step:(cv + dop_range)]; %#ok<AGROW>
        end
        all_dops = unique(round(all_dops / dop_step) * dop_step);

        all_dels = [];
        for ci = 1:length(coarse_ests.cand_delays)
            cd = coarse_ests.cand_delays(ci);
            all_dels = [all_dels, (cd-1):(cd+1)]; %#ok<AGROW>
        end
        all_dels = unique(all_dels);
        all_dels = all_dels(all_dels >= 0 & all_dels <= derivedParams.padLen);  % delays beyond the fixed guard are infeasible

        [gd, gv] = ndgrid(all_dels, all_dops);
        dict_d = gd(:)';  dict_v = gv(:)';
        n_atoms = length(dict_d);

        dict_ests = struct('AoAs', zeros(1,n_atoms), 'AoDs', zeros(1,n_atoms), ...
                           'pathDelays', dict_d, 'pathDopplers', dict_v);
        Phi_dict = build_response_matrix(dict_ests, params, derivedParams);

        if params.verbose_output
            fprintf('  SSR dict (CLEAN-seeded): %d delays × %d dopplers = %d atoms\n', ...
                length(all_dels), length(all_dops), n_atoms);
        end

        selected = block_omp(Phi_dict, y, blk_rows, n_blocks, K);

        init_ests = struct('AoAs', zeros(1,K), 'AoDs', zeros(1,K), ...
                           'pathDelays', dict_d(selected), ...
                           'pathDopplers', dict_v(selected));

        if params.verbose_output
            fprintf('  SSR block-OMP: delay=[%s] Doppler=[%s]\n', ...
                sprintf('%.0f ', init_ests.pathDelays), ...
                sprintf('%.2f ', init_ests.pathDopplers));
        end

    elseif nargin >= 6 && ~isempty(coarse_ests) && isfield(coarse_ests, 'pathDelays')
        %% === Paired coarse estimates (legacy) ===
        init_ests = coarse_ests;
    else
        %% === No coarse: full-grid OMP ===
        all_delays   = 0:1:derivedParams.padLen;
        all_dopplers = 0:1:(N-1);

        n_atoms = length(all_delays) * length(all_dopplers);
        [grid_r, grid_v] = ndgrid(all_delays, all_dopplers);
        cand_delays   = grid_r(:).';
        cand_dopplers = grid_v(:).';
        cand_ests = struct('AoAs', zeros(1,n_atoms), 'AoDs', zeros(1,n_atoms), ...
                           'pathDelays', cand_delays, 'pathDopplers', cand_dopplers);
        Phi_dict = build_response_matrix(cand_ests, params, derivedParams);

        if params.verbose_output
            fprintf('  SSR dict: %d atoms (full grid)\n', n_atoms);
        end

        selected = block_omp(Phi_dict, y, blk_rows, n_blocks, K);

        init_ests = struct('AoAs', zeros(1,K), 'AoDs', zeros(1,K), ...
                           'pathDelays', round(cand_delays(selected)), ...
                           'pathDopplers', cand_dopplers(selected));
    end

    %% === Finalize: wrap Dopplers + per-antenna LS ===
    final_ests = init_ests;
    final_ests.pathDopplers = mod(final_ests.pathDopplers + params.N/2, params.N) - params.N/2;
    final_ests.beta_per_tx = per_antenna_ls(YTF, XTF, final_ests, params, derivedParams, n0);

    if params.verbose_output
        fprintf('--- SSR final ---\n');
        fprintf('SSR. Delay:   %s\n', sprintf('%5.2f ', final_ests.pathDelays));
        fprintf('SSR. Doppler: %s\n', sprintf('%5.2f ', final_ests.pathDopplers));
        fprintf('---------------------------------\n');
    end
end


%% ========================================================================
%  Block-OMP: per-antenna correlation + per-antenna LS subtraction
%  ========================================================================

function selected = block_omp(Phi, y, blk_rows, n_blocks, K)
%BLOCK_OMP OMP with per-block (per-antenna) coefficients.
%   Correlation: sum_b |Phi_b' * y_resid_b|^2
%   Subtraction: per-block LS with all selected atoms.

    n_atoms = size(Phi, 2);

    Phi_blk = cell(n_blocks, 1);
    y_blk   = cell(n_blocks, 1);
    y_resid = cell(n_blocks, 1);
    for b = 1:n_blocks
        Phi_blk{b} = Phi(blk_rows{b}, :);
        y_blk{b}   = y(blk_rows{b});
        y_resid{b} = y_blk{b};
    end

    selected = zeros(1, K);
    for k = 1:K
        corr_spec = zeros(n_atoms, 1);
        for b = 1:n_blocks
            corr_spec = corr_spec + abs(Phi_blk{b}' * y_resid{b}).^2;
        end

        if k > 1
            corr_spec(selected(1:k-1)) = 0;
        end

        [~, best] = max(corr_spec);
        selected(k) = best;

        for b = 1:n_blocks
            Phi_sel = Phi_blk{b}(:, selected(1:k));
            alpha_b = Phi_sel \ y_blk{b};
            y_resid{b} = y_blk{b} - Phi_sel * alpha_b;
        end
    end
end


%% ========================================================================
%  Per-antenna LS: recover per-TX channel gains
%  ========================================================================

function beta_per_tx = per_antenna_ls(YTF, XTF, ests, params, derivedParams, n0)
%PER_ANTENNA_LS For each TX antenna, solve ridge LS using its pilot blocks.
%   Pilot model y = c*Phi*beta + noise, ridge s2 (pilot_noise_model).

    K  = length(ests.pathDelays);
    Nt = params.Nt;
    pilot_blocks = params.pilot_blocks;

    est0 = struct('AoAs', zeros(1,K), 'AoDs', zeros(1,K), ...
                  'pathDelays', ests.pathDelays, ...
                  'pathDopplers', ests.pathDopplers);

    [c, s2] = pilot_noise_model(XTF, params, n0);
    beta_per_tx = zeros(K, Nt);
    for tx = 1:Nt
        blks_tx = pilot_blocks([pilot_blocks.tx] == tx);
        if isempty(blks_tx), continue; end

        y_tx = build_measurement_vector(YTF, XTF, blks_tx, params.Nr, n0);
        params_tx = params;
        params_tx.pilot_blocks = blks_tx;
        Phi_tx = c * build_response_matrix(est0, params_tx, derivedParams);

        beta_per_tx(:, tx) = (Phi_tx'*Phi_tx + s2*eye(K)) \ (Phi_tx'*y_tx);
    end
end


%% ========================================================================
%  NMSE: per-TX DD channel reconstruction
%  ========================================================================

function nmse_val = compute_nmse_nr1(chanParams, ests, params, derivedParams)
    HDD_true = generate_DD_channel(chanParams, params, derivedParams);

    HDD_est = generate_DD_channel(per_tx_channel_nr1(ests), params, derivedParams);

    nmse_val = sum(abs(HDD_true - HDD_est).^2, 'all') / sum(abs(HDD_true).^2, 'all');
end


function est = per_tx_channel_nr1(ests)
% Estimated channel with per-Tx gains (paths-by-Nt) in place of the Tx steering.
    est = struct('AoAs', zeros(1, length(ests.pathDelays)), ...
                 'pathDelays', ests.pathDelays, 'pathDopplers', ests.pathDopplers, ...
                 'pathGains', ests.beta_per_tx, 'gain_per_tx', true);
end


%% ========================================================================
%  Per-target error reference (errors: shared compute_errors, gated by ssr.range)
%  ========================================================================

function truth = cluster_means_nr1(chanParams)
% Cluster-mean parameters, the reference for per-target error statistics.
    truth = struct('AoAs', zeros(1, chanParams.n_clusters), ...
        'AoDs', chanParams.mean_AoDs, 'pathDelays', chanParams.mean_delays, ...
        'pathDopplers', chanParams.mean_dopplers);
end


%% ========================================================================
%  3GPP SCM clustered channel generation
%  ========================================================================

function chanParams = generate_channel_params_nr1(params)
%GENERATE_CHANNEL_PARAMS_NR1 3GPP SCM clustered channel (TR 25.996).
%   Ns=1: point target. Ns=20: 3GPP with fixed AoD offsets per cluster.

    K  = params.targets_num;
    Ns = 1;
    angular_spread = 2;
    if isfield(params, 'Ns_subpaths'),       Ns = params.Ns_subpaths; end
    if isfield(params, 'angular_spread_deg'), angular_spread = params.angular_spread_deg; end

    alpha_table = [0.0447, 0.1413, 0.2492, 0.3715, 0.5129, ...
                   0.6797, 0.8844, 1.1481, 1.5195, 2.1551];

    base = generate_channel_params(params);  % cluster delays, Dopplers and gains
    mean_delays   = base.pathDelays;
    mean_dopplers = base.pathDopplers;

    mean_AoDs = params.target_aods(1:K);

    total = K * Ns;
    chanParams.pathDelays   = zeros(1, total);
    chanParams.pathDopplers = zeros(1, total);
    chanParams.AoDs         = zeros(1, total);
    chanParams.AoAs         = zeros(1, total);
    chanParams.pathGains    = zeros(1, total);

    for i = 1:K
        idx = (i-1)*Ns + (1:Ns);

        chanParams.pathDelays(idx)   = mean_delays(i);
        chanParams.pathDopplers(idx) = mean_dopplers(i);

        % 3GPP TR 25.996 Table 5.2: fixed subpath AoD offsets (20 = 10 pairs of +/-)
        if Ns == 20
            offsets_deg = angular_spread * ...
                reshape([alpha_table; -alpha_table], 1, []);
            chanParams.AoDs(idx) = mean_AoDs(i) + offsets_deg;
        elseif Ns > 1
            n_pairs = min(ceil(Ns/2), 10);
            off = angular_spread * ...
                reshape([alpha_table(1:n_pairs); -alpha_table(1:n_pairs)], 1, []);
            chanParams.AoDs(idx) = mean_AoDs(i) + off(1:Ns);
        else
            chanParams.AoDs(idx) = mean_AoDs(i);
        end

        % Equal-power subpaths with random phases, total cluster power = 1
        g_cluster = base.pathGains(i);
        chanParams.pathGains(idx) = (g_cluster / sqrt(Ns)) * exp(1j * 2*pi * rand(1, Ns));
    end

    chanParams.n_clusters    = K;
    chanParams.Ns_subpaths   = Ns;
    chanParams.mean_AoDs     = mean_AoDs;
    chanParams.mean_delays   = mean_delays;
    chanParams.mean_dopplers = mean_dopplers;

    if params.verbose_output
        fprintf('--- Channel (Nt=%d, Nr=1, Ns=%d) ---\n', params.Nt, Ns);
        fprintf('True Delay:    %s\n', sprintf('%5.2f ', mean_delays));
        fprintf('True Doppler:  %s\n', sprintf('%5.2f ', mean_dopplers));
        fprintf('True AoD:      %s\n', sprintf('%5.2f ', mean_AoDs));
        if Ns > 1
            fprintf('Angular spread: %.1f deg, Total subpaths: %d\n', angular_spread, total);
        end
        fprintf('---------------------------------\n');
    end
end
