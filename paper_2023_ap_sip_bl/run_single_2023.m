%RUN_SINGLE_2023 Run, summarize, plot, and optionally save AP-SIP CSI results.

close all;
params = config_params_2023();
results = run_sim_2023(params);

mean_nmse_pa = mean(results.NMSE_PA, 1);
mean_ber_pa = mean(results.BER_PA, 1);
compute_ber = results.params.compute_ber;
fprintf('\n--- 2023 AP-SIP PA-BL Summary ---\n');
for snr_idx = 1:numel(results.SNR_db_list)
    if compute_ber
        fprintf('SNR %g dB: PA-BL operator NMSE %.2f dB, PA-BL BER %.4e\n', ...
            results.SNR_db_list(snr_idx), ...
            10 * log10(mean_nmse_pa(snr_idx)), mean_ber_pa(snr_idx));
    else
        fprintf('SNR %g dB: PA-BL operator NMSE %.2f dB, BER skipped.\n', ...
            results.SNR_db_list(snr_idx), 10 * log10(mean_nmse_pa(snr_idx)));
    end
end

if params.trials > 1
    out_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    out_file = fullfile(out_dir, build_result_filename(params));
    save(out_file, 'results');
    fprintf('Saved: %s\n', out_file);
end

figure;
if compute_ber
    subplot(1, 2, 1);
end
plot(results.SNR_db_list, 10 * log10(mean_nmse_pa), '-o', ...
    'LineWidth', 2, 'DisplayName', 'PA-BL operator NMSE');
grid on;
xlabel('SNR [dB]');
ylabel('Channel NMSE [dB]');
title('PA-BL CSI Estimation');
legend('Location', 'best');

if compute_ber
    subplot(1, 2, 2);
    semilogy(results.SNR_db_list, mean_ber_pa, '-s', 'LineWidth', 2);
    grid on;
    xlabel('SNR [dB]');
    ylabel('PA-BL BER');
    title('Data Detection');
end
