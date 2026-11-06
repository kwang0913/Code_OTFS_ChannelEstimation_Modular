function results = run_ber_batch(overrides, config_fn)
%RUN_BER_BATCH Multi-trial BER run with data detection on every trial.
%   RESULTS = RUN_BER_BATCH(OVERRIDES) runs run_sim_v2 on
%   config_params_std(OVERRIDES) with compute_ber = true and saves the
%   results to results/ber_<build_result_filename>.
%   RESULTS = RUN_BER_BATCH(OVERRIDES, CONFIG_FN) uses CONFIG_FN instead,
%   e.g. @config_params_HF.
%
%   Bash: matlab -batch "run_ber_batch(struct('trials', 50))"

    if nargin < 1 || isempty(overrides), overrides = struct(); end
    if nargin < 2, config_fn = @config_params_std; end
    overrides.compute_ber = true;
    params = config_fn(overrides);
    results = run_sim_v2(params);

    out_dir = 'results';
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    out_file = fullfile(out_dir, ['ber_' build_result_filename(params)]);
    save(out_file, 'results');
    fprintf('Saved: %s\n', out_file);
end
