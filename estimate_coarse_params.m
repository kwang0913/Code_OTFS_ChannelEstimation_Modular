function coarse_ests = estimate_coarse_params(YTF, XTF, chanParams, params, n0, derivedParams)
%ESTIMATE_COARSE_PARAMS CLEAN AoA + per-angle DFT for delay/Doppler + law-of-cosines AoD.

    interpolate_len = 64;

    pilot_blocks = params.pilot_blocks;
    delay_blks   = pilot_blocks(strcmp({pilot_blocks.dim}, 'delay'));
    doppler_blks = pilot_blocks(strcmp({pilot_blocks.dim}, 'doppler'));

    %% 1-2. Coarse AoA via CLEAN successive cancellation
    rx_OTFS_clean = reshape(YTF, [], params.Nr);
    x_ang = asind(((0:interpolate_len-1) - interpolate_len/2) / (params.ddlr * interpolate_len));
    max_targets = params.targets_num;
    orig_power = norm(rx_OTFS_clean, 'fro')^2;

    detected_angles = [];
    original_peak = 0;
    Cor_angle_init = [];

    for iter = 1:max_targets
        % Spatial FFT with interpolation
        rx_idft = fft(rx_OTFS_clean, interpolate_len, 2);
        Cor_angle = median(abs(rx_idft), 1);
        Cor_angle = fftshift(Cor_angle);

        if iter == 1
            original_peak = max(Cor_angle);
            Cor_angle_init = Cor_angle / original_peak;
        end

        spec = Cor_angle / original_peak;

        % Find prominent peaks
        [pks, locs, ~, proms] = findpeaks(Cor_angle, 'SortStr', 'descend');

        % Stop if no peaks, or prominence < 10% of original (residual is noise-like)
        if isempty(pks)
            break;
        end
        if proms(1) / original_peak < 0.1 && iter > 1
            break;
        end

        peak_val = pks(1);
        peak_idx = locs(1);

        % Stop if spectrum is flat: peak barely above median, left-skewed noise floor
        spec_mean = mean(spec);
        spec_med  = median(spec);
        if spec_mean < spec_med && (peak_val/original_peak - spec_med) < (spec_med - spec_mean) && iter > 1
            break;
        end

        theta_k = x_ang(peak_idx);
        detected_angles(end+1) = theta_k; %#ok<AGROW>

        % Steering vector for detected angle
        a_k = exp(1j * 2 * pi * (0:params.Nr-1)' * params.ddlr * sind(theta_k));

        % LS estimate per sample and subtract
        alpha_k = (a_k' * a_k) \ (a_k' * rx_OTFS_clean.');
        rx_OTFS_clean = rx_OTFS_clean - (a_k * alpha_k).';
    end

    n_detected = length(detected_angles);

    % Fallback: if nothing detected, use strongest bin
    if n_detected == 0
        rx_idft = fft(reshape(YTF, [], params.Nr), interpolate_len, 2);
        Cor_angle = median(abs(rx_idft), 1);
        Cor_angle = fftshift(Cor_angle);
        [~, fb_idx] = max(Cor_angle);
        detected_angles = x_ang(fb_idx);
        n_detected = 1;
        Cor_angle_init = Cor_angle / max(Cor_angle);
    end

    ang = detected_angles;

    if params.show_figures
        figure; hold on;
        h = plot(x_ang, Cor_angle_init, '-', 'DisplayName', 'Coarse AoAs', 'LineWidth', 2);
        for j_idx = 1:params.targets_num
            xline(chanParams.AoAs(j_idx), 'r--', 'LineWidth', 2);
        end
        hold off;
        legend(h, 'Location', 'southwest', 'FontSize', 16);
        xlabel('Angle [\circ]', 'FontSize', 16); ylabel('Normalized Magnitude', 'FontSize', 16);
        set(gca, 'FontSize', 14, 'FontWeight', 'Bold');

        % Zoomed inset
        inset_ax = axes('Position', [0.25, 0.6, 0.3, 0.25]);
        hold(inset_ax, 'on');
        plot(inset_ax, x_ang, Cor_angle_init, '-', 'LineWidth', 2);
        for j_idx = 1:params.targets_num
            xline(inset_ax, chanParams.AoAs(j_idx), 'r--', 'LineWidth', 2);
        end
        hold(inset_ax, 'off');
        [~, peak_idx] = max(Cor_angle_init);
        peak_ang = x_ang(peak_idx);
        xlim(inset_ax, [peak_ang-5, peak_ang+5]);
        ylim(inset_ax, [0.8, 1.1]);
        set(inset_ax, 'XTick', [], 'YTick', []);
        box(inset_ax, 'on'); grid(inset_ax, 'on');
        set(inset_ax, 'FontSize', 12);
    end

    % Each CLEAN detection is one target (no grouping needed)
    ang_profile = cell(1, n_detected);
    for i = 1:n_detected
        ang_profile{i} = i;
    end

    Cor.AoA = ang;

    %% 3. Per-angle-bin amplitude profile (MMSE beamforming)
    Ang_mtx = exp(1j * 2 * pi * (0:params.Nr-1)' * params.ddlr * sind(round(ang,1)));
    W_bf = (Ang_mtx' * Ang_mtx + n0 * eye(length(ang))) \ Ang_mtx';
    ATF = zeros(length(ang), params.M, params.N);
    for m_idx = 1:params.M
        for n_idx = 1:params.N
            ATF(:, m_idx, n_idx) = W_bf * squeeze(YTF(m_idx, n_idx, :));
        end
    end

    %% 4. Coarse delay — coherent combining across delay blocks
    interpolate_len = 128;
    x_delay = (0:interpolate_len-1) * params.M / interpolate_len;
    Cor.delay = zeros(1, n_detected);

    % Coherent combining: place pilot ratios at their subcarrier offsets so
    % a single IFFT sees the full aperture. (n, tx) factors cancel in the ratio.
    max_offset_delay = 0;
    for db = 1:length(delay_blks)
        max_offset_delay = max(max_offset_delay, max(delay_blks(db).offsets));
    end
    n_combined_delay = max_offset_delay + 1;

    r_combined = zeros(n_combined_delay, length(ang));
    for db = 1:length(delay_blks)
        blk = delay_blks(db);
        for p = 1:blk.n_pilots
            for j_idx = 1:length(ang)
                r_combined(blk.offsets(p) + 1, j_idx) = ...
                    (ATF(j_idx, blk.m(p), blk.n) / XTF(blk.m(p), blk.n, blk.tx)) / ...
                    (ATF(j_idx, blk.ref_m, blk.ref_n) / XTF(blk.ref_m, blk.ref_n, blk.ref_tx));
            end
        end
    end

    r_idft_combined = abs(ifft(r_combined .* hann(n_combined_delay), interpolate_len, 1));

    for p_idx = 1:length(ang)
        if isempty(ang_profile{p_idx}), continue; end

        temp = r_idft_combined(:, p_idx);
        Cor_delay_profile = temp / max(temp);

        if params.show_figures
            % Find nearest true target by angle for overlay
            [~, nearest_true] = min(abs(chanParams.AoAs - ang(p_idx)));
            true_delay_val = chanParams.pathDelays(nearest_true);

            figure; hold on;
            h_plot = plot(x_delay, Cor_delay_profile, '-', 'LineWidth', 2);
            xline(true_delay_val, 'r--', 'LineWidth', 2);
            legend(h_plot, 'Coarse Delay', 'Location', 'southwest', 'FontSize', 16);
            xlabel('Delay [tap]', 'FontSize', 16); ylabel('Normalized Magnitude', 'FontSize', 16);
            set(gca, 'FontSize', 14, 'FontWeight', 'Bold'); hold off;

            % Zoomed inset
            inset_ax = axes('Position', [0.55, 0.6, 0.3, 0.25]);
            hold(inset_ax, 'on');
            plot(inset_ax, x_delay, Cor_delay_profile, '-', 'LineWidth', 2);
            xline(inset_ax, true_delay_val, 'r--', 'LineWidth', 2);
            hold(inset_ax, 'off');
            xlim(inset_ax, [true_delay_val-6, true_delay_val+6]);
            ylim(inset_ax, [0.8, 1.1]);
            set(inset_ax, 'XTick', [], 'YTick', []);
            box(inset_ax, 'on'); grid(inset_ax, 'on');
            set(inset_ax, 'FontSize', 12);
        end

        n_expected = length(ang_profile{p_idx});
        [~, idx] = findpeaks(Cor_delay_profile, 'MinPeakProminence', 0.01, ...
            'SortStr', 'descend', 'NPeaks', n_expected);
        if isempty(idx), [~, idx] = max(Cor_delay_profile); end
        idx = sort(idx);

        for j_idx = 1:length(ang_profile{p_idx})
            if j_idx <= length(idx)
                Cor.delay(ang_profile{p_idx}(j_idx)) = x_delay(idx(j_idx));
            else
                Cor.delay(ang_profile{p_idx}(j_idx)) = x_delay(idx(end));
            end
        end
    end

    %% 5. Coarse Doppler — coherent combining across doppler blocks
    x_doppler = -((0:interpolate_len-1) - interpolate_len/2) * params.N / interpolate_len;
    Cor.doppler = zeros(1, n_detected);

    % Coherent combining: place pilot ratios at symbol offsets for Doppler IFFT.
    max_offset_doppler = 0;
    for db = 1:length(doppler_blks)
        max_offset_doppler = max(max_offset_doppler, max(doppler_blks(db).offsets));
    end
    n_combined_doppler = max_offset_doppler + 1;

    v_combined = zeros(n_combined_doppler, length(ang));
    for db = 1:length(doppler_blks)
        blk = doppler_blks(db);
        for p = 1:blk.n_pilots
            for j_idx = 1:length(ang)
                v_combined(blk.offsets(p) + 1, j_idx) = ...
                    (ATF(j_idx, blk.m, blk.n(p)) / XTF(blk.m, blk.n(p), blk.tx)) / ...
                    (ATF(j_idx, blk.ref_m, blk.ref_n) / XTF(blk.ref_m, blk.ref_n, blk.ref_tx));
            end
        end
    end

    v_idft_combined = abs(ifft(v_combined .* hann(n_combined_doppler), interpolate_len, 1));

    for p_idx = 1:length(ang)
        if isempty(ang_profile{p_idx}), continue; end

        temp = v_idft_combined(:, p_idx);
        temp_shifted = [temp(interpolate_len/2+1:end); temp(1:interpolate_len/2)];
        Cor_doppler_profile = temp_shifted / max(temp_shifted);

        if params.show_figures
            % Find nearest true target by angle for overlay
            [~, nearest_true] = min(abs(chanParams.AoAs - ang(p_idx)));
            true_doppler_val = chanParams.pathDopplers(nearest_true);

            figure; hold on;
            h_plot = plot(x_doppler, Cor_doppler_profile, '-', 'LineWidth', 2);
            xline(true_doppler_val, 'r--', 'LineWidth', 2);
            legend(h_plot, 'Coarse Doppler', 'Location', 'southwest', 'FontSize', 16);
            xlabel('Doppler [tap]', 'FontSize', 16); ylabel('Normalized Magnitude', 'FontSize', 16);
            set(gca, 'FontSize', 14, 'FontWeight', 'Bold'); hold off;

            % Zoomed inset
            inset_ax = axes('Position', [0.55, 0.6, 0.3, 0.25]);
            hold(inset_ax, 'on');
            plot(inset_ax, x_doppler, Cor_doppler_profile, '-', 'LineWidth', 2);
            xline(inset_ax, true_doppler_val, 'r--', 'LineWidth', 2);
            hold(inset_ax, 'off');
            xlim(inset_ax, [true_doppler_val-3, true_doppler_val+3]);
            ylim(inset_ax, [0.8, 1.1]);
            set(inset_ax, 'XTick', [], 'YTick', []);
            box(inset_ax, 'on'); grid(inset_ax, 'on');
            set(inset_ax, 'FontSize', 12);
        end

        n_expected = length(ang_profile{p_idx});
        [~, idx] = findpeaks(Cor_doppler_profile, 'MinPeakProminence', 0.01, ...
            'SortStr', 'descend', 'NPeaks', n_expected);
        if isempty(idx), [~, idx] = max(Cor_doppler_profile); end
        idx = sort(idx);

        for j_idx = 1:length(ang_profile{p_idx})
            if j_idx <= length(idx)
                Cor.doppler(ang_profile{p_idx}(j_idx)) = x_doppler(idx(j_idx));
            else
                Cor.doppler(ang_profile{p_idx}(j_idx)) = x_doppler(idx(end));
            end
        end
    end

    %% 6. Coarse AoD via Law of Cosines
    range = Cor.delay / (params.M * params.df) * params.c;
    R0 = params.init.R0;
    Cor.RoA = (range.^2 - R0^2) ./ (2 * (range - R0 * cosd(abs(Cor.AoA))));
    Cor.RoD = range - Cor.RoA;

    cos_val = (R0^2 + Cor.RoD.^2 - Cor.RoA.^2) ./ (2 * R0 * Cor.RoD);
    cos_val(cos_val > 1) = 1; cos_val(cos_val < -1) = -1;
    Cor.AoD = acosd(cos_val) .* sign(Cor.AoA);

    if params.show_figures
        figure; hold on; grid on;
        stem(Cor.AoD, ones(1, length(Cor.AoD)), '-', 'LineWidth', 2);
        h = 'Coarse AoDs';
        for j_idx = 1:params.targets_num
            xline(chanParams.AoDs(j_idx), 'r--', 'LineWidth', 2);
        end
        hold off;
        legend(h, 'Location', 'southwest', 'FontSize', 16);
        xlabel('Angle [\circ]', 'FontSize', 16); ylabel('Normalized Magnitude', 'FontSize', 16);
        set(gca, 'FontSize', 14, 'FontWeight', 'Bold');
    end

    %% 7. Package results
    coarse_ests.AoAs = round(Cor.AoA, 1);
    coarse_ests.pathDelays = round(Cor.delay);
    coarse_ests.pathDopplers = round(Cor.doppler,1);
    coarse_ests.AoDs = round(Cor.AoD, 1);
    coarse_ests.n_detected = n_detected;

    %% 8. Estimate path gains via LS
    Phi_coarse = build_response_matrix(coarse_ests, params);
    K = n_detected;
    Nr = params.Nr;
    n_meas = compute_n_meas(params.pilot_blocks, Nr);
    y = zeros(n_meas, 1);
    idx = 0;
    for b = 1:length(params.pilot_blocks)
        blk = params.pilot_blocks(b);
        y(idx + (1:Nr)) = squeeze(YTF(blk.ref_m, blk.ref_n, 1:Nr)) ./ XTF(blk.ref_m, blk.ref_n, blk.ref_tx);
        idx = idx + Nr;
        for p = 1:blk.n_pilots
            if strcmp(blk.dim, 'delay')
                y(idx + (1:Nr)) = squeeze(YTF(blk.m(p), blk.n, 1:Nr)) ./ XTF(blk.m(p), blk.n, blk.tx);
            else
                y(idx + (1:Nr)) = squeeze(YTF(blk.m, blk.n(p), 1:Nr)) ./ XTF(blk.m, blk.n(p), blk.tx);
            end
            idx = idx + Nr;
        end
    end
    coarse_ests.pathGains = (Phi_coarse'*Phi_coarse + n0*eye(K)) \ (Phi_coarse'*y);

    if params.verbose_output
        fprintf('--- Coarse Estimation Results ---\n');
        fprintf('Cor. AoA:     %s\n', sprintf('%5.2f ', coarse_ests.AoAs));
        fprintf('Cor. AoD:     %s\n', sprintf('%5.2f ', coarse_ests.AoDs));
        fprintf('Cor. Delay:   %s\n', sprintf('%5.2f ', coarse_ests.pathDelays));
        fprintf('Cor. Doppler: %s\n', sprintf('%5.2f ', coarse_ests.pathDopplers));
        fprintf('---------------------------------\n');
    end
end
