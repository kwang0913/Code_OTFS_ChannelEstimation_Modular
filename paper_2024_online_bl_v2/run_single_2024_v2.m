%RUN_SINGLE_2024_V2 Run, summarize, plot, and save main-physical ORGBL results.
close all;
params = config_params_2024_v2();
results = run_sim_2024_v2(params);
mean_nmse = mean(results.NMSE_ORGBL, 1);
skipped = all(results.detection_skipped(:));
if ~skipped, mean_ber = mean(results.BER_ORGBL, 1); end
fprintf('\n--- 2024 v2 Main-Physical ORGBL Summary ---\n');
for index = 1:numel(results.SNR_db_list)
    if skipped
        fprintf('SNR %g dB: main ORGBL NMSE %.2f dB, BER skipped.\n', ...
            results.SNR_db_list(index), 10*log10(mean_nmse(index)));
    else
        fprintf('SNR %g dB: main ORGBL NMSE %.2f dB, BER %.4e\n', ...
            results.SNR_db_list(index), 10*log10(mean_nmse(index)), ...
            mean_ber(index));
    end
end
if params.trials > 1
    output_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
    if exist(output_dir, 'dir') ~= 7, mkdir(output_dir); end
    filename = sprintf('online_bl_2024_v2_M%d_N%d_Nt%d_Nr%d_trials%d.mat', ...
        params.M, params.N, params.Nt, params.Nr, params.trials);
    save(fullfile(output_dir, filename), 'results', '-v7.3');
end
figure;
if skipped
    plot(results.SNR_db_list, 10*log10(mean_nmse), '-o', 'LineWidth', 2);
    grid on; xlabel('SNR [dB]'); ylabel('Main physical NMSE [dB]');
    title('2024 v2 Main-Physical ORGBL');
else
    subplot(1,2,1);
    plot(results.SNR_db_list, 10*log10(mean_nmse), '-o', 'LineWidth', 2);
    grid on; xlabel('SNR [dB]'); ylabel('Main physical NMSE [dB]');
    title('Channel Estimation');
    subplot(1,2,2);
    semilogy(results.SNR_db_list, mean_ber, '-s', 'LineWidth', 2);
    grid on; xlabel('SNR [dB]'); ylabel('BER'); title('Data Detection');
end
