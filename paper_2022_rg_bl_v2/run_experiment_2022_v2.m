function run_experiment_2022_v2(json_file)
%RUN_EXPERIMENT_2022_V2 Run JSON-configured Srivastava 2022 sweeps.
    if ~(ischar(json_file) || (isstring(json_file) && isscalar(json_file)))
        error('run_experiment_2022_v2:InvalidFile','json_file must be text.');
    end
    json_file = char(json_file);
    if exist(json_file,'file') ~= 2
        error('run_experiment_2022_v2:MissingFile','Missing JSON file: %s',json_file);
    end
    config = jsondecode(fileread(json_file));
    if isfield(config,'overrides'), base = config.overrides; else, base = struct(); end
    if isfield(config,'sweep'), sweep = config.sweep; else, sweep = []; end
    run_count = max(1,numel(sweep));
    folder = fileparts(mfilename('fullpath'));
    addpath(fileparts(folder),folder);
    output_dir = fullfile(folder,'results');
    if exist(output_dir,'dir') ~= 7, mkdir(output_dir); end
    [~,experiment_name] = fileparts(json_file);
    for run_index = 1:run_count
        overrides = base;
        if ~isempty(sweep)
            if iscell(sweep), item = sweep{run_index}; else, item = sweep(run_index); end
            names = fieldnames(item);
            for index = 1:numel(names), overrides.(names{index}) = item.(names{index}); end
        end
        params = config_params_2022_v2(overrides);
        results = run_sim_2022_v2(params);
        if run_count == 1
            filename = experiment_name+".mat";
        else
            filename = sprintf('%s_run%d.mat',experiment_name,run_index);
        end
        save(fullfile(output_dir,filename),'results','-v7.3');
        fprintf('Saved: %s\n',fullfile(output_dir,filename));
    end
end
