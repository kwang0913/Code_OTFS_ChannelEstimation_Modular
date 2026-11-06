function results = run_ber_batch_2023(overrides)
%RUN_BER_BATCH_2023 Multi-trial 2023 run with data detection on every trial.
%   RESULTS = RUN_BER_BATCH_2023(OVERRIDES) runs run_sim_2023 on
%   config_params_2023(OVERRIDES) with compute_ber = true and saves the
%   results to paper_2023_ap_sip_bl/results/ber_<build_result_filename>.
    if nargin < 1 || isempty(overrides), overrides = struct(); end
    overrides.compute_ber = true;
    params = config_params_2023(overrides);
    results = run_sim_2023(params);

    out_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    out_file = fullfile(out_dir, ['ber_' build_result_filename(params)]);
    save(out_file, 'results');
    fprintf('Saved: %s\n', out_file);
end
