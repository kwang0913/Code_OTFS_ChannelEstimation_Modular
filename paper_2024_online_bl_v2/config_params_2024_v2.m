function params = config_params_2024_v2(overrides)
%CONFIG_PARAMS_2024_V2 Configure the main-physical 2024 online BL model.
    if nargin < 1 || isempty(overrides)
        overrides = struct();
    end
    if ~isstruct(overrides) || ~isscalar(overrides)
        error('config_params_2024_v2:InvalidOverrides', ...
            'overrides must be a scalar struct.');
    end

    folder = fileparts(mfilename('fullpath'));
    project_root = fileparts(folder);
    baseline_folder = fullfile(project_root, 'paper_2024_online_bl');
    addpath(project_root, '-begin');
    addpath(baseline_folder, '-begin');
    params = config_params_2024(overrides);

    max_delay = max(params.delay_grid_2024(:));
    cp_length = max_delay;
    if isfield(overrides, 'cp_length_2024')
        cp_length = overrides.cp_length_2024;
    end
    if ~isscalar(cp_length) || ~isreal(cp_length) || ~isfinite(cp_length) || ...
            cp_length < max_delay || cp_length ~= round(cp_length)
        error('config_params_2024_v2:InvalidCPLength', ...
            'cp_length_2024 must be an integer covering the full delay grid.');
    end

    params.padType = 'CP';
    params.cp_length_2024 = cp_length;
    params.padding_offset_2024 = 0;
    params.doppler_reference_columns_2024 = params.N;
    params.total_columns_2024 = params.N + params.N_pilots_2024;
end
