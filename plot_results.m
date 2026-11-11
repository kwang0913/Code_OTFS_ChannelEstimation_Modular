function plot_results(results)
%PLOT_RESULTS NMSE, BER, conditional MSE and detection tables from run_sim results.
%   plot_results(results) accepts the struct returned by run_sim or loaded
%   from a saved .mat file. Failed trials stay NaN and are ignored.

    params = results.params;
    snr = results.SNR_db_list;
    K = params.targets_num;
    ok = ~results.failed;

    % --- NMSE (and BER for single-trial runs) ---
    figure;
    subplot(1, 2, 1);
    plot(snr, 10*log10(mean(results.nmse_coarse, 1, 'omitnan')), '-o', 'LineWidth', 2); hold on;
    plot(snr, 10*log10(mean(results.nmse_ssr, 1, 'omitnan')), '-s', 'LineWidth', 2);
    grid on; xlabel('SNR [dB]'); ylabel('NMSE [dB]'); legend('Coarse', 'SSR');
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');
    if params.trials == 1
        subplot(1, 2, 2);
        ber = results.BER(1, :);
        zero_ber = ber == 0;  % measured zero, labelled below
        ber(zero_ber) = NaN;  % a logarithmic axis cannot represent zero
        semilogy(snr, ber, '-s', 'LineWidth', 2);
        grid on; xlabel('SNR [dB]'); ylabel('BER'); legend('SSR + LMMSE');
        set(gca, 'FontSize', 12, 'FontWeight', 'bold');
        if any(zero_ber)
            if all(zero_ber), ylim([1e-6, 1]); end
            text(0.03, 0.07, sprintf('BER = 0 at SNR [dB]: %s\n(omitted on log axis)', ...
                strtrim(sprintf('%g ', snr(zero_ber)))), 'Units', 'normalized', ...
                'VerticalAlignment', 'bottom', 'FontSize', 9);
        end
    end

    % --- MSE (conditional on detection) ---
    labels = {'Angle [deg]', 'Delay [taps]', 'Doppler'};
    figure;
    for p = 1:3
        subplot(1, 3, p);
        mse_coarse = squeeze(mean(mean(results.err_coarse(:, :, :, p).^2, 1, 'omitnan'), 3, 'omitnan'));
        mse_ssr    = squeeze(mean(mean(results.err_ssr(:, :, :, p).^2, 1, 'omitnan'), 3, 'omitnan'));
        semilogy(snr, mse_coarse, '-o', 'LineWidth', 2); hold on;
        semilogy(snr, mse_ssr, '-s', 'LineWidth', 2);
        grid on; xlabel('SNR [dB]'); ylabel(['MSE | matched - ' labels{p}]); legend('Coarse', 'SSR');
        set(gca, 'FontSize', 11, 'FontWeight', 'bold');
    end
    sgtitle('Conditional MSE (matched reflections only)');

    % --- Detection: reflections matched to a true target within ssr.range ---
    fprintf('\n--- Detection (%d targets, matched within ssr.range) ---\n', K);
    fprintf('%-16s', 'SNR(dB)');
    fprintf('%8d ', snr);
    fprintf('\n');
    rows = {'n_detected (coarse)', results.n_detected_coarse; 'n_detected (SSR)', results.n_detected_ssr; ...
            'P_d (coarse)', results.n_matched_coarse / K; 'P_d (SSR)', results.n_matched_ssr / K};
    for i = 1:size(rows, 1)
        v = rows{i, 2};  v(~ok) = NaN;
        fprintf('%-20s', rows{i, 1});
        fprintf('%8.2f ', mean(v, 1, 'omitnan'));
        fprintf('\n');
    end
end
