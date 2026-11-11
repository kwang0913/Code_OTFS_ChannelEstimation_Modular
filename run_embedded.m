% run_embedded.m
% Embedded DD pilot estimation benchmark. Pilots in DD domain with guard zones.

close all;
clear all;

params = config_params_std(struct('selected_targets', [1,2,3,4]));
derivedParams = compute_derived_params(params);
use_parallel = params.trials > 1 && params.n_workers > 0;

NMSE = zeros(params.trials, length(params.SNR_db_list));

% Guard sizes from channel spread (delay: fixed padLen, which covers the largest total delay)
max_delay = params.padLen;
max_doppler = ceil(max(abs([params.init.pathDopplers, params.doppler_range])));

% Pilot locations in DD domain (one per TX, spaced at M/8 intervals to avoid overlap)
pilot_locations = [
    1*params.M/8, params.N/2;
    3*params.M/8, params.N/2;
    5*params.M/8, params.N/2;
    7*params.M/8, params.N/2
];

for snr_idx = 1:length(params.SNR_db_list)
    SNR_db = params.SNR_db_list(snr_idx);
    fprintf('Running SNR = %d dB\n', SNR_db);

    temp_nmse = zeros(params.trials, 1);
    dq = parallel.pool.DataQueue;
    afterEach(dq, @(~) update_progress());

    if use_parallel
        params.show_figures = false;
        params.verbose_output = false;
        parfor (trial_idx = 1:params.trials, params.n_workers)
            rng(params.rng_seed_start + trial_idx - 1, 'threefry');
            temp_nmse(trial_idx) = run_embedded_trial(params, derivedParams, pilot_locations, ...
                max_delay, max_doppler, SNR_db);
            send(dq, 0);
        end
    else
        for trial_idx = 1:params.trials
            rng(params.rng_seed_start + trial_idx - 1, 'threefry');
            temp_nmse(trial_idx) = run_embedded_trial(params, derivedParams, pilot_locations, ...
                max_delay, max_doppler, SNR_db);
            send(dq, 0);
        end
    end
    fprintf('\n');
    NMSE(:,snr_idx) = temp_nmse;
    fprintf('SNR %d dB done.\n', SNR_db);
end

% --- Save Results ---
results.NMSE           = NMSE;
results.SNR_db_list    = params.SNR_db_list;
results.params         = params;
results.pilot_locations = pilot_locations;

out_dir = 'results';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
[~, base_noext] = fileparts(build_result_filename(params));
out_file = fullfile(out_dir, ['embed_' base_noext '.mat']);
save(out_file, 'results');
fprintf('Saved: %s\n', out_file);

figure;
plot(params.SNR_db_list, 10*log10(mean(NMSE, 1)), '-o', 'LineWidth', 2);
grid on;
xlabel('SNR [dB]');
ylabel('NMSE [dB]');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

function nmse_val = run_embedded_trial(params, derivedParams, pilot_locations, max_delay, max_doppler, SNR_db)
    chanParams = generate_channel_params(params);
    HDD_true = generate_DD_channel(chanParams, params, derivedParams);

    % Random QAM data grid
    databits = randi([0 1], params.M*log2(params.ModOrder), params.N, params.Nt);
    datagrid = zeros(params.M, params.N, params.Nt);
    for tx_idx = 1:params.Nt
        temp = bit2int(squeeze(databits(:,:,tx_idx)), log2(params.ModOrder));
        if isfield(params, 'modulation') && strcmpi(params.modulation, 'PSK')
            datagrid(:,:,tx_idx) = pskmod(temp, params.ModOrder);
        else
            datagrid(:,:,tx_idx) = qammod(temp, params.ModOrder, UnitAveragePower=true);
        end
    end

    % Zero guard zone around each pilot (delay +/-1, Doppler +/-2 beyond spread), then place pilot
    XDD_pilot = datagrid;
    for tx_idx = 1:params.Nt
        loc = pilot_locations(tx_idx, :);
        XDD_pilot(loc(1)-max_delay-1 : loc(1)+max_delay+1, ...
                  loc(2)-2*max_doppler-2 : loc(2)+2*max_doppler+2, :) = 0;
        XDD_pilot(loc(1), loc(2), tx_idx) = params.pilot_power;
    end

    % DD -> TF, pass through channel
    XTF_pilot = ifft(fft(XDD_pilot, [], 1), [], 2) / sqrt(params.M/params.N);
    [~, YDD_received, ~] = channel_and_noise_rect(XTF_pilot, chanParams, params, derivedParams, SNR_db);

    if params.show_figures
        figure;
        mesh(0:params.N-1, 0:params.M-1, abs(YDD_received(:,:,1)));
        xlabel('Doppler'); ylabel('Delay'); zlabel('|Y_{DD}|');
        title(sprintf('Received DD Domain (Rx 1, SNR = %d dB)', SNR_db));
        view([-30 60]);
        set(gca, 'FontSize', 12, 'FontWeight', 'bold');
    end

    % CFAR-like threshold: estimate noise sigma from bins outside guard zones
    guard_mask = false(params.M, params.N);
    for tx_idx = 1:params.Nt
        loc = pilot_locations(tx_idx, :);
        guard_mask(loc(1)-max_delay-1 : loc(1)+max_delay+1, ...
                   loc(2)-2*max_doppler-2 : loc(2)+2*max_doppler+2) = true;
    end
    outside = abs(YDD_received);
    outside(repmat(guard_mask, 1, 1, params.Nr)) = 0;
    outside_vals = outside(outside > 0);
    sigma_hat = sqrt(mean(outside_vals.^2) / 2);  % Rayleigh MLE for noise std
    P_fa = 1e-4;
    tap_threshold = sigma_hat * sqrt(-2 * log(P_fa));  % Neyman-Pearson threshold

    % Extract channel estimate from received DD around each pilot
    HDD_est = zeros(params.M, params.N, params.Nr, params.Nt);
    for tx_idx = 1:params.Nt
        loc = pilot_locations(tx_idx, :);
        delay_indices   = loc(1)-max_delay : loc(1)+max_delay;
        doppler_indices = loc(2)-2*max_doppler-1 : loc(2)+2*max_doppler+1;

        for rx_idx = 1:params.Nr
            temp = zeros(params.M, params.N);
            temp(delay_indices, doppler_indices) = YDD_received(delay_indices, doppler_indices, rx_idx);
            [lp, vp] = find(abs(temp) > tap_threshold);

            % Expand each detected tap to include Doppler neighbours [vp-5, vp+5]
            if ~isempty(lp)
                lp_exp = [];
                vp_exp = [];
                for kk = 1:length(lp)
                    offsets = -1:1;
                    vp_neighbours = vp(kk) + offsets;
                    % Clip to valid range [1, N]
                    vp_neighbours = vp_neighbours(vp_neighbours >= 1 & vp_neighbours <= params.N);
                    lp_exp = [lp_exp; repmat(lp(kk), length(vp_neighbours), 1)];
                    vp_exp = [vp_exp; vp_neighbours(:)];
                end
                % Merge with originals and remove duplicates
                pairs = unique([lp_exp, vp_exp], 'rows');
                lp = pairs(:,1);
                vp = pairs(:,2);
            end

            if params.show_figures
                peaks_only = zeros(params.M, params.N);
                peaks_only(sub2ind(size(peaks_only), lp, vp)) = abs(temp(sub2ind(size(temp), lp, vp)));
                figure;
                mesh(0:params.N-1, 0:params.M-1, peaks_only);
                xlabel('Doppler'); ylabel('Delay'); zlabel('|temp|');
                title(sprintf('Detected Taps (Tx %d, Rx %d, SNR = %d dB)', tx_idx, rx_idx, SNR_db));
                view([-30 60]);
                set(gca, 'FontSize', 12, 'FontWeight', 'bold');
            end

            % Map each detected tap back to DD channel estimate (undo pilot offset, mod wrap)
            for k = 1:length(lp)
                tap_val = temp(lp(k), vp(k)) / params.pilot_power;
                delay_offset = lp(k) - loc(1);
                doppler_offset = vp(k) - loc(2);
                dest_delay = mod(delay_offset, params.M) + 1;
                dest_doppler = mod(doppler_offset, params.N) + 1;
                HDD_est(dest_delay, dest_doppler, rx_idx, tx_idx) = tap_val;
            end
        end
    end

    numerator = sum(abs(HDD_true - HDD_est).^2, 'all');
    denominator = sum(abs(HDD_true).^2, 'all');
    nmse_val = numerator / denominator;
end

function update_progress()
    persistent count t0
    if isempty(count)
        count = 0;
        t0 = tic;
    end
    count = count + 1;
    fprintf('\r  Progress: %d (%.1f s)', count, toc(t0));
end
