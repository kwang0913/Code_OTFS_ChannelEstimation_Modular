%RUN_SINGLE_2023 Run, summarize, plot, and optionally save AP-SIP CSI results.

close all;
params = config_params_2023();
results = run_sim_2023(params);

mean_nmse_pa = mean(results.NMSE_PA, 1);
mean_ber_pa = mean(results.BER_PA, 1);
fprintf('\n--- 2023 AP-SIP PA-BL Summary ---\n');
for snr_idx = 1:numel(results.SNR_db_list)
    fprintf('SNR %g dB: PA-BL operator NMSE %.2f dB, PA-BL BER %.4e\n', ...
        results.SNR_db_list(snr_idx), ...
        10 * log10(mean_nmse_pa(snr_idx)), mean_ber_pa(snr_idx));
end

if params.trials > 1
    out_dir = 'results';
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    out_name = sprintf('ap_sip_bl_2023_M%d_N%d_Nt%d_Nr%d_trials%d.mat', ...
        params.M, params.N, params.Nt, params.Nr, params.trials);
    out_file = fullfile(out_dir, out_name);
    save(out_file, 'results');
    fprintf('Saved: %s\n', out_file);
end

figure;
subplot(1, 2, 1);
plot(results.SNR_db_list, 10 * log10(mean_nmse_pa), '-o', ...
    'LineWidth', 2, 'DisplayName', 'PA-BL operator NMSE');
grid on;
xlabel('SNR [dB]');
ylabel('Channel NMSE [dB]');
title('PA-BL CSI Estimation');
legend('Location', 'best');

subplot(1, 2, 2);
semilogy(results.SNR_db_list, mean_ber_pa, '-s', 'LineWidth', 2);
grid on;
xlabel('SNR [dB]');
ylabel('PA-BL BER');
title('Data Detection');

