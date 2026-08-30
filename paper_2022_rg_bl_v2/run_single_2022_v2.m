%RUN_SINGLE_2022_V2 Run, summarize, plot, and optionally save 2022 results.
close all;
clearvars;
params = config_params_2022_v2();
results = run_sim_2022_v2(params);
mean_rgbl = mean(results.NMSE_RGBL,1);
mean_rgomp = mean(results.NMSE_RGOMP,1);
mean_ber_rgbl = mean(results.BER_RGBL,1);
mean_ber_rgomp = mean(results.BER_RGOMP,1);
fprintf('\n--- Srivastava 2022 Main-Physical Summary ---\n');
for index = 1:numel(results.SNR_db_list)
    fprintf(['SNR %g dB: RG-BL NMSE %.2f dB, RG-OMP NMSE %.2f dB, ' ...
        'RG-BL BER %.4e, RG-OMP BER %.4e\n'],results.SNR_db_list(index), ...
        10*log10(mean_rgbl(index)),10*log10(mean_rgomp(index)), ...
        mean_ber_rgbl(index),mean_ber_rgomp(index));
end
if params.trials ~= 1
    output_dir = fullfile(fileparts(mfilename('fullpath')),'results');
    if exist(output_dir,'dir') ~= 7, mkdir(output_dir); end
    filename = sprintf('rg_bl_2022_v2_M%d_N%d_Nt%d_Nr%d_trials%d.mat', ...
        params.M,params.N,params.Nt,params.Nr,params.trials);
    save(fullfile(output_dir,filename),'results','-v7.3');
end
figure;
subplot(1,2,1);
plot(results.SNR_db_list,10*log10(mean_rgbl),'-o','LineWidth',2, ...
    'DisplayName','RG-BL'); hold on;
plot(results.SNR_db_list,10*log10(mean_rgomp),'--s','LineWidth',2, ...
    'DisplayName','RG-OMP'); hold off; grid on; legend('Location','best');
xlabel('SNR [dB]'); ylabel('Main physical NMSE [dB]'); title('Channel Estimation');
subplot(1,2,2);
semilogy(results.SNR_db_list,mean_ber_rgbl,'-o','LineWidth',2, ...
    'DisplayName','RG-BL'); hold on;
semilogy(results.SNR_db_list,mean_ber_rgomp,'--s','LineWidth',2, ...
    'DisplayName','RG-OMP'); hold off; grid on; legend('Location','best');
xlabel('SNR [dB]'); ylabel('BER'); title('Data Detection');
