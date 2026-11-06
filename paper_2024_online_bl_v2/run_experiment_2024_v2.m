function run_experiment_2024_v2(json_file)
%RUN_EXPERIMENT_2024_V2 Run JSON-configured main-physical 2024 experiments.
    if ~(ischar(json_file) || (isstring(json_file) && isscalar(json_file)))
        error('run_experiment_2024_v2:InvalidFile', ...
            'json_file must be a character vector or string scalar.');
    end
    json_file = char(json_file);
    if exist(json_file, 'file') ~= 2
        error('run_experiment_2024_v2:MissingFile', ...
            'JSON experiment file does not exist: %s', json_file);
    end
    config = jsondecode(fileread(json_file));
    if ~isstruct(config) || ~isscalar(config)
        error('run_experiment_2024_v2:InvalidConfig', ...
            'The JSON root must decode to a scalar object.');
    end
    base = struct();
    if isfield(config, 'overrides'), base = config.overrides; end
    if ~isstruct(base) || ~isscalar(base)
        error('run_experiment_2024_v2:InvalidConfig', ...
            'overrides must be a JSON object.');
    end
    sweep = {};
    if isfield(config, 'sweep'), sweep = config.sweep; end
    run_count = max(1, numel(sweep));
    [~, experiment_name] = fileparts(json_file);
    folder = fileparts(mfilename('fullpath'));
    addpath(fileparts(folder), folder);
    output_dir = fullfile(folder, 'results');
    if exist(output_dir, 'dir') ~= 7, mkdir(output_dir); end
    for run_index = 1:run_count
        overrides = base;
        if ~isempty(sweep)
            if iscell(sweep), entry = sweep{run_index}; else, entry = sweep(run_index); end
            overrides = apply_overrides(overrides, entry);
        end
        params = config_params_2024_v2(overrides);
        results = run_sim_2024_v2(params);
        if run_count == 1
            filename = sprintf('%s.mat', experiment_name);
        else
            filename = sprintf('%s_run%d.mat', experiment_name, run_index);
        end
        output_file = fullfile(output_dir, filename);
        save(output_file, 'results', '-v7.3');
        fprintf('Saved: %s\n', output_file);
    end
end
