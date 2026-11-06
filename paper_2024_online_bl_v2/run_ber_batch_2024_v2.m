function results = run_ber_batch_2024_v2(overrides)
%RUN_BER_BATCH_2024_V2 Multi-trial 2024 run with data detection on every trial.
%   RESULTS = RUN_BER_BATCH_2024_V2(OVERRIDES) runs run_sim_2024_v2 on
%   config_params_2024_v2(OVERRIDES) with compute_ber = true and saves the
%   results to paper_2024_online_bl_v2/results/ber_<build_result_filename>.
%   Saved with -v7.3: per-trial records exceed the 2 GB v7 limit at scale.
    if nargin < 1 || isempty(overrides), overrides = struct(); end
    overrides.compute_ber = true;
    params = config_params_2024_v2(overrides);
    results = run_sim_2024_v2(params);

    out_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    out_file = fullfile(out_dir, ['ber_' build_result_filename(params)]);
    save(out_file, 'results', '-v7.3');
    fprintf('Saved: %s\n', out_file);
end
