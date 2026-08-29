% run_single.m
% Main entry script. Runs one simulation, displays results and plots.

close all;
clear all;

% --- Simulation Control ---
% params = config_params_HF();
params = config_params_std();

% --- Run Simulation ---
% results = run_sim(params);
results = run_sim_v2(params);

% --- Save Results (only for multi-trial runs) ---
if params.trials > 1
    out_dir = 'results';
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    out_file = fullfile(out_dir, build_result_filename(params));
    save(out_file, 'results');
    fprintf('Saved: %s\n', out_file);
end

% --- Plot Final Results ---
figure;
subplot(1,2,1);
plot(results.SNR_db_list, 10*log10(mean(results.NMSE_before, 1)), '-o', 'LineWidth', 2); hold on;
plot(results.SNR_db_list, 10*log10(mean(results.NMSE_after, 1)), '-s', 'LineWidth', 2);
grid on;
xlabel('SNR [dB]');
ylabel('NMSE [dB]');
legend('Coarse', 'Final (SSR)');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

% --- Display Coarse AoA Detection ---
fprintf('\n--- Coarse AoA Detection ---\n');
fprintf('%-12s', 'SNR(dB)');
fprintf('%8d ', results.SNR_db_list);
fprintf('\n');
fprintf('%-12s', 'n_detected');
fprintf('%8.2f ', mean(results.n_detected_aoa, 1));
fprintf('\n');
fprintf('%-12s', 'P_d');
fprintf('%8.2f ', mean(results.n_detected_aoa == params.targets_num, 1));
fprintf('\n');

% % --- Display Per-Target Errors (mean +/- std across trials) ---
% param_labels = {'AoA', 'AoD', 'Delay', 'Doppler', 'Gain'};
% K = size(results.pct_err_coarse, 3);
% nSNR = length(results.SNR_db_list);
% fprintf('\n--- Per-Target Abs Errors (mean +/- std across trials) ---\n');
% fprintf('%-16s', 'SNR(dB)');
% fprintf('%16d ', results.SNR_db_list);
% fprintf('\n');
% fprintf('  Coarse:\n');
% for k = 1:K
%     fprintf('    Target %d:\n', k);
%     for p = 1:5
%         fprintf('    %-12s', param_labels{p});
%         for s = 1:nSNR
%             m = mean(results.pct_err_coarse(:, s, k, p), 1, 'omitnan');
%             sd = std(results.pct_err_coarse(:, s, k, p), 0, 1, 'omitnan');
%             fprintf('%7.2f+/-%-7.2f', m, sd);
%         end
%         fprintf('\n');
%     end
% end
% fprintf('  SSR:\n');
% for k = 1:K
%     fprintf('    Target %d:\n', k);
%     for p = 1:5
%         fprintf('    %-12s', param_labels{p});
%         for s = 1:nSNR
%             m = mean(results.pct_err_ssr(:, s, k, p), 1, 'omitnan');
%             sd = std(results.pct_err_ssr(:, s, k, p), 0, 1, 'omitnan');
%             fprintf('%7.2f+/-%-7.2f', m, sd);
%         end
%         fprintf('\n');
%     end
% end
% fprintf('--------------------------------------\n\n');

% --- Display Per-Target Coarse Estimates (mean +/- std across trials) ---
param_labels = {'AoA', 'AoD', 'Delay', 'Doppler', 'Gain'};
K = size(results.coarse_ests_all, 3);
nSNR = length(results.SNR_db_list);
fprintf('--- Per-Target Coarse Estimates (mean +/- std across trials) ---\n');
fprintf('%-16s', 'SNR(dB)');
fprintf('%16d ', results.SNR_db_list);
fprintf('\n');
for k = 1:K
    fprintf('  Target %d (true: DoA=%.1f, DoD=%.1f, delay=%.1f, doppler=%.1f):\n', ...
        k, params.init.DoAs(k), params.init.DoDs(k), ...
        params.init.pathDelays(k), params.init.pathDopplers(k));
    for p = 1:5
        fprintf('    %-12s', param_labels{p});
        for s = 1:nSNR
            m = mean(results.coarse_ests_all(:, s, k, p), 1, 'omitnan');
            sd = std(results.coarse_ests_all(:, s, k, p), 0, 1, 'omitnan');
            fprintf('%7.2f+/-%-7.2f', m, sd);
        end
        fprintf('\n');
    end
end
fprintf('--------------------------------------\n\n');

if any(results.BER_after(:) > 0)
    subplot(1,2,2);
    semilogy(results.SNR_db_list, mean(results.BER_after, 1), '-s', 'LineWidth', 2);
    grid on;
    xlabel('SNR [dB]');
    ylabel('BER');
    legend('Final (SSR)');
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');
end
