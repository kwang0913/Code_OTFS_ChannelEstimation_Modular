function coarse_ests = estimate_coarse_params_v2(YTF, XTF, chanParams, params, n0, derivedParams)
%ESTIMATE_COARSE_PARAMS_V2 Simplified coarse: CLEAN AoA, beamform, per-angle DFT, geometry AoD.
%   1. AoA CLEAN  2. MMSE beamform  3. Per-angle delay/Doppler DFT  4. AoD via law-of-cosines

    pilot_blocks = params.pilot_blocks;
    delay_blks   = pilot_blocks(strcmp({pilot_blocks.dim}, 'delay'));
    doppler_blks = pilot_blocks(strcmp({pilot_blocks.dim}, 'doppler'));

    %% 1. Setup
    delay_blk = delay_blks(1);
    b_d = delay_blk.n_pilots;
    stride_d = 1;
    if b_d > 1, stride_d = delay_blk.m(2) - delay_blk.m(1); end
    P_delay = 2048;  % oversampled DFT for fine-grained delay peak detection
    delay_axis = (0:P_delay-1) * params.M / (P_delay * stride_d);
    valid_d = delay_axis < params.M / stride_d;

    doppler_blk = doppler_blks(1);
    b_v = doppler_blk.n_pilots;
    stride_v = 1;
    if b_v > 1, stride_v = doppler_blk.n(2) - doppler_blk.n(1); end
    P_doppler = 2048;
    doppler_axis = ((0:P_doppler-1) - P_doppler/2) * params.N / (P_doppler * stride_v);

    x_delay = XTF(delay_blk.m, delay_blk.n, delay_blk.tx);
    x_doppler = squeeze(XTF(doppler_blk.m, doppler_blk.n, doppler_blk.tx));

    %% 2. AoA CLEAN
    interpolate_len = 64;
    max_targets = params.targets_num;
    rx_OTFS_clean = reshape(YTF, [], params.Nr);
    x_ang = asind(((0:interpolate_len-1) - interpolate_len/2) / (params.ddlr * interpolate_len));

    detected_angles = [];
    original_peak = 0;
    Cor_angle_init = [];

    for iter = 1:max_targets
        rx_idft = fft(rx_OTFS_clean, interpolate_len, 2);
        Cor_angle = fftshift(median(abs(rx_idft), 1));

        if iter == 1
            original_peak = max(Cor_angle);
            Cor_angle_init = Cor_angle / original_peak;
        end

        [pks, locs, ~, proms] = findpeaks(Cor_angle, 'SortStr', 'descend');
        if isempty(pks), break; end
        if proms(1) / original_peak < 0.1 && iter > 1, break; end

        spec = Cor_angle / original_peak;
        spec_mean = mean(spec);  spec_med = median(spec);
        if spec_mean < spec_med && (pks(1)/original_peak - spec_med) < (spec_med - spec_mean) && iter > 1
            break;
        end

        theta_k = x_ang(locs(1));
        detected_angles(end+1) = theta_k; %#ok<AGROW>

        a_k = exp(1j * 2 * pi * (0:params.Nr-1)' * params.ddlr * sind(theta_k));
        alpha_k = (a_k' * a_k) \ (a_k' * rx_OTFS_clean.');
        rx_OTFS_clean = rx_OTFS_clean - (a_k * alpha_k).';
    end

    n_detected = length(detected_angles);
    if n_detected == 0
        rx_idft = fft(reshape(YTF, [], params.Nr), interpolate_len, 2);
        Cor_angle = fftshift(median(abs(rx_idft), 1));
        [~, fb_idx] = max(Cor_angle);
        detected_angles = x_ang(fb_idx);
        n_detected = 1;
        Cor_angle_init = Cor_angle / max(Cor_angle);
    end
    ang = detected_angles;

    if params.show_figures
        figure; hold on;
        plot(x_ang, Cor_angle_init, '-', 'DisplayName', 'Coarse AoAs', 'LineWidth', 2);
        for j = 1:params.targets_num
            xline(chanParams.AoAs(j), 'r--', 'LineWidth', 2);
        end
        hold off; legend('Location', 'southwest', 'FontSize', 16);
        xlabel('Angle [\circ]', 'FontSize', 16); ylabel('Normalized Magnitude', 'FontSize', 16);
        set(gca, 'FontSize', 14, 'FontWeight', 'Bold');
    end

    %% 3. Beamform + per-angle DFT
    Ang_mtx = exp(1j * 2 * pi * (0:params.Nr-1)' * params.ddlr * sind(round(ang, 1)));
    W_bf = (Ang_mtx' * Ang_mtx + n0 * eye(n_detected)) \ Ang_mtx';
    % Reshape YTF to (Nr x M*N), apply beamformer, reshape back to (n_detected x M x N)
    ATF = reshape(W_bf * reshape(permute(YTF, [3,1,2]), params.Nr, []), n_detected, params.M, params.N);

    Cor.AoA     = ang;
    Cor.delay   = zeros(1, n_detected);
    Cor.doppler = zeros(1, n_detected);

    for p = 1:n_detected
        % Delay: IFFT peak
        ratio_d = squeeze(ATF(p, delay_blk.m, delay_blk.n)).' ./ x_delay;
        spec_d = abs(ifft(ratio_d, P_delay));
        spec_d(~valid_d') = 0;
        [~, pk] = max(spec_d);
        Cor.delay(p) = delay_axis(pk);

        % Doppler: FFT peak
        ratio_v = squeeze(ATF(p, doppler_blk.m, doppler_blk.n)) ./ x_doppler(:);
        spec_v = abs(fftshift(fft(ratio_v, P_doppler)));
        [~, pk] = max(spec_v);
        Cor.doppler(p) = doppler_axis(pk);
    end
    Cor.delay = round(Cor.delay);

    if params.show_figures
        figure;
        for p = 1:min(n_detected, 4)
            ratio_d = squeeze(ATF(p, delay_blk.m, delay_blk.n)).' ./ x_delay;
            spec = abs(ifft(ratio_d, P_delay)); spec(~valid_d') = 0; spec = spec/max(spec);
            subplot(2,2,p); hold on;
            plot(delay_axis(valid_d), spec(valid_d), '-', 'LineWidth', 1.5);
            [~, nt] = min(abs(chanParams.AoAs - ang(p)));
            xline(chanParams.pathDelays(nt), 'r--', 'LineWidth', 1.5);
            xline(Cor.delay(p), 'g:', 'LineWidth', 2);
            hold off; grid on; title(sprintf('%.1f°', ang(p)));
            xlabel('Delay [tap]'); set(gca, 'FontSize', 10, 'FontWeight', 'bold');
        end
        sgtitle('Per-Angle Delay DFT', 'FontWeight', 'bold');

        figure;
        for p = 1:min(n_detected, 4)
            ratio_v = squeeze(ATF(p, doppler_blk.m, doppler_blk.n)) ./ x_doppler(:);
            spec = abs(fftshift(fft(ratio_v, P_doppler))); spec = spec/max(spec);
            subplot(2,2,p); hold on;
            plot(doppler_axis, spec, '-', 'LineWidth', 1.5);
            [~, nt] = min(abs(chanParams.AoAs - ang(p)));
            xline(chanParams.pathDopplers(nt), 'r--', 'LineWidth', 1.5);
            xline(Cor.doppler(p), 'g:', 'LineWidth', 2);
            hold off; grid on; title(sprintf('%.1f°', ang(p)));
            xlabel('Doppler [tap]'); set(gca, 'FontSize', 10, 'FontWeight', 'bold');
        end
        sgtitle('Per-Angle Doppler DFT', 'FontWeight', 'bold');
    end

    %% 4. AoD + path gains
    [Cor.AoD, ~] = aod_from_geometry(Cor.AoA, Cor.delay, params.M, params.df, params.c, params.init.R0);

    coarse_ests.AoAs         = round(Cor.AoA, 1);
    coarse_ests.pathDelays   = round(Cor.delay);
    coarse_ests.pathDopplers = round(Cor.doppler, 1);
    coarse_ests.AoDs         = round(Cor.AoD, 1);
    coarse_ests.n_detected   = n_detected;

    Phi_coarse = build_response_matrix(coarse_ests, params);
    y = build_measurement_vector(YTF, XTF, params.pilot_blocks, params.Nr);
    coarse_ests.pathGains = (Phi_coarse'*Phi_coarse + n0*eye(n_detected)) \ (Phi_coarse'*y);

    if params.verbose_output
        fprintf('--- Coarse Estimation Results ---\n');
        fprintf('Cor. AoA:     %s\n', sprintf('%5.2f ', coarse_ests.AoAs));
        fprintf('Cor. AoD:     %s\n', sprintf('%5.2f ', coarse_ests.AoDs));
        fprintf('Cor. Delay:   %s\n', sprintf('%5.2f ', coarse_ests.pathDelays));
        fprintf('Cor. Doppler: %s\n', sprintf('%5.2f ', coarse_ests.pathDopplers));
        fprintf('---------------------------------\n');
    end
end
