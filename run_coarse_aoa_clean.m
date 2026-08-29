% run_coarse_aoa_clean.m
% Standalone CLEAN AoA test. Visualizes successive cancellation iterations.
% Finds the strongest spatial peak, LS-estimates and subtracts its
% contribution, then repeats on the residual.

close all; clear;

%% Setup
params = config_params_std(struct('Nr', 16, 'show_figures', false));

SNR_db = 20;

[chanParams, derivedParams] = generate_channel_params(params);
[XTF, ~, params] = generate_TF(params);
[YTF, ~, ~] = channel_and_noise_rect(XTF, chanParams, params, derivedParams, SNR_db);

fprintf('True AoAs: %s\n', sprintf('%.1f° ', chanParams.AoAs));

%% CLEAN parameters
interpolate_len = 128;
max_targets = params.targets_num+1;
n0 = 1 / 10^(SNR_db / 10);

%% Angle axis (same convention as estimate_coarse_params)
x_ang = asind(((0:interpolate_len-1) - interpolate_len/2) / (params.ddlr * interpolate_len));

%% Successive cancellation
rx_OTFS = reshape(YTF, [], params.Nr);   % (M*N) x Nr
orig_power = norm(rx_OTFS, 'fro')^2;

detected_angles = [];
detected_peaks  = [];
spectra = {};
original_peak = 0;

for iter = 1:max_targets
    % Spatial FFT with interpolation
    rx_idft   = fft(rx_OTFS, interpolate_len, 2);
    Cor_angle = median(abs(rx_idft), 1);
    Cor_angle = fftshift(Cor_angle);

    if iter == 1
        original_peak = max(Cor_angle);
    end

    % Store normalized spectrum for this iteration
    spec = Cor_angle / original_peak;
    spectra{iter} = spec;

    % Residual power ratio
    res_ratio = norm(rx_OTFS, 'fro')^2 / orig_power;

    % Find prominent peaks
    [pks, locs, ~, proms] = findpeaks(Cor_angle, 'SortStr', 'descend');
    % proms is available for debugging; not used in stopping criteria

    % Stopping criterion 1: prominence < 10% of original
    if isempty(pks)
        fprintf('Iter %d: no peaks found — stopping\n', iter);
        break;
    end
    if proms(1) / original_peak < 0.1 && iter > 1
        fprintf('Iter %d: prom/original = %.4f < 0.1 — stopping (prominence)\n', ...
            iter, proms(1) / original_peak);
        break;
    end

    peak_val = pks(1);
    peak_idx = locs(1);

    % Stopping criterion 2: spectrum is flat (mean < median, peak doesn't stand out)
    spec_mean = mean(spec);
    spec_med  = median(spec);
    if spec_mean < spec_med && (peak_val/original_peak - spec_med) < (spec_med - spec_mean) && iter > 1
        fprintf('Iter %d: spectrum flat (mean=%.4f < med=%.4f) — stopping\n', ...
            iter, spec_mean, spec_med);
        break;
    end

    theta_k = x_ang(peak_idx);
    detected_angles(end+1) = theta_k;
    detected_peaks(end+1)  = peak_val / original_peak;

    fprintf('Iter %d: angle = %6.2f°, rel. peak = %.4f, prom = %.4f, residual = %.4f\n', ...
        iter, theta_k, detected_peaks(end), proms(1), res_ratio);

    % Steering vector for detected angle
    a_k = exp(1j * 2 * pi * (0:params.Nr-1)' * params.ddlr * sind(theta_k));

    % LS estimate per sample and subtract
    alpha_k = (a_k' * a_k) \ (a_k' * rx_OTFS.');
    rx_OTFS = rx_OTFS - (a_k * alpha_k).';
end

%% Plot — one subplot per iteration
n_det  = length(detected_angles);
n_cols = min(n_det, 4);
n_rows = ceil(n_det / n_cols);

figure('Position', [100 100 400*n_cols 350*n_rows]);
for iter = 1:n_det
    subplot(n_rows, n_cols, iter);
    hold on;

    % Residual spectrum
    plot(x_ang, spectra{iter}, 'b-', 'LineWidth', 1.5);

    % Mark detected peak
    [~, pidx] = min(abs(x_ang - detected_angles(iter)));
    plot(x_ang(pidx), spectra{iter}(pidx), 'ro', 'MarkerSize', 10, 'LineWidth', 2);

    % Mean and median spectrum level
    yline(mean(spectra{iter}), 'g--', 'LineWidth', 1.5);
    yline(median(spectra{iter}), 'm--', 'LineWidth', 1.5);

    % True AoAs
    for j = 1:params.targets_num
        xline(chanParams.AoAs(j), 'r--', 'LineWidth', 1.5);
    end

    hold off;
    title(sprintf('Iter %d: \\theta = %.1f°', iter, detected_angles(iter)));
    xlabel('Angle [°]'); ylabel('Norm. Magnitude');
    grid on; ylim([0 1.1]);
    set(gca, 'FontSize', 11);
end
sgtitle(sprintf('CLEAN AoA (SNR = %d dB, Nr = %d)', SNR_db, params.Nr), ...
    'FontSize', 14, 'FontWeight', 'Bold');

%% Summary
fprintf('\n--- Summary ---\n');
fprintf('Detected: %s\n', sprintf('%6.1f° ', detected_angles));
fprintf('True AoA: %s\n', sprintf('%6.1f° ', chanParams.AoAs));
