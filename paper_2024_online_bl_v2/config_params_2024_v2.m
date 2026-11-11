function params = config_params_2024_v2(overrides)
%CONFIG_PARAMS_2024_V2 Configure the main-physical 2024 online BL model.
    if nargin < 1 || isempty(overrides), overrides = struct(); end
    if ~isstruct(overrides) || ~isscalar(overrides)
        error('config_params_2024_v2:InvalidOverrides', ...
            'overrides must be a scalar struct.');
    end
    folder = fileparts(mfilename('fullpath'));
    project_root = fileparts(folder);
    addpath(project_root, '-begin');
    % Online pilots are independent of the root delay-Doppler pilot layout; the
    % root only sizes a random-pilot layout, which keeps small grids valid.
    base_overrides = overrides;
    if ~isfield(overrides, 'pilot_blocks_def'), base_overrides.pilot_blocks_def = []; end
    params = config_params_std(base_overrides);
    params.config_name = 'online_bl_2024_v2';

    defaults = struct( ...
        'SNR_db_list', -10:5:30, 'trials', 2, 'use_parallel', false, ...
        'compute_ber', [], ...
        'verbose_output', true, ...
        'N_pilots_2024', 4, ...
        'doppler_oversampling_2024', 4, 'data_power_2024', 1, ...
        'pilot_power_2024', 1, 'paper_gain_mode', 'fading', ...
        'em_tol_2024', 1e-6, 'em_max_iter_2024', 50, ...
        'variance_floor_2024', 1e-12, ...
        'lmmse_tol_2024', 1e-8, 'lmmse_max_iter_2024', 200);
    names = fieldnames(defaults);
    for index = 1:numel(names)
        if ~isfield(overrides, names{index})
            params.(names{index}) = defaults.(names{index});
        end
    end
    [configured_delays, configured_dopplers] = local_configured_targets(params);
    required_delay_count = ceil(max(configured_delays)) + 1;
    if ~isfield(overrides, 'M_tau'), params.M_tau = required_delay_count; end
    local_positive_integer(params.M_tau, 'M_tau');
    if params.M_tau < required_delay_count || params.M_tau > params.M
        error('config_params_2024_v2:DelayBounds', ...
            'M_tau must cover configured delays and not exceed M.');
    end
    params.delay_grid_2024 = (0:params.M_tau-1).';

    oversampling = local_positive_integer( ...
        params.doppler_oversampling_2024, 'doppler_oversampling_2024');
    if isempty(params.doppler_range)
        nu_max = max(abs(configured_dopplers));
    else
        if numel(params.doppler_range) ~= 2 || ...
                any(~isfinite(params.doppler_range(:))) || ...
                any(params.doppler_range(:) < 0) || ...
                params.doppler_range(1) > params.doppler_range(2)
            error('config_params_2024_v2:InvalidDopplerRange', ...
                'doppler_range must be an ascending nonnegative two-element vector.');
        end
        nu_max = max(abs(params.doppler_range));
    end
    required_doppler_count = 2*ceil(nu_max*oversampling)+1;
    if ~isfield(overrides, 'G_nu'), params.G_nu = required_doppler_count; end
    local_positive_integer(params.G_nu, 'G_nu');
    if mod(params.G_nu, 2) ~= 1 || params.G_nu < required_doppler_count
        error('config_params_2024_v2:DopplerBounds', ...
            'G_nu must be odd and cover the configured symmetric Doppler range.');
    end
    center = (params.G_nu-1)/2;
    params.doppler_grid_2024 = ((0:params.G_nu-1)-center).' / oversampling;

    integer_names = {'N_pilots_2024','em_max_iter_2024','lmmse_max_iter_2024'};
    for index = 1:numel(integer_names)
        local_positive_integer(params.(integer_names{index}), integer_names{index});
    end
    if params.N_pilots_2024 > params.N
        error('config_params_2024_v2:PilotCount', ...
            'N_pilots_2024 must not exceed N.');
    end
    local_positive_scalar(params.data_power_2024, 'data_power_2024');
    local_positive_scalar(params.pilot_power_2024, 'pilot_power_2024');
    local_positive_scalar(params.variance_floor_2024, 'variance_floor_2024');
    local_nonnegative_scalar(params.em_tol_2024, 'em_tol_2024');
    local_positive_scalar(params.lmmse_tol_2024, 'lmmse_tol_2024');

    max_delay = max(params.delay_grid_2024);
    cp_length = max_delay;
    if isfield(overrides, 'cp_length_2024'), cp_length = overrides.cp_length_2024; end
    if ~isscalar(cp_length) || ~isreal(cp_length) || ~isfinite(cp_length) || ...
            cp_length < max_delay || cp_length > params.M || ...
            cp_length ~= round(cp_length)
        error('config_params_2024_v2:InvalidCPLength', ...
            'cp_length_2024 must cover the delay grid and not exceed M.');
    end
    params.padType = 'CP';
    params.cp_length_2024 = cp_length;
    params.padding_offset_2024 = 0;
    params.doppler_reference_columns_2024 = params.N;
    params.total_columns_2024 = params.N + params.N_pilots_2024;
end

function [delays, dopplers] = local_configured_targets(params)
    if ~isempty(params.n_targets) && isfield(params, 'target_pool')
        source = params.target_pool;
        delays = [source.DoD] + [source.DoA];
        dopplers = [source.Doppler];
    else
        delays = params.init.DoDs + params.init.DoAs;
        dopplers = params.init.pathDopplers;
    end
    if isempty(delays) || isempty(dopplers) || any(~isfinite(delays)) || ...
            any(delays < 0) || any(~isfinite(dopplers))
        error('config_params_2024_v2:InvalidTargets', ...
            'Configured delays and Dopplers must be finite and nonempty.');
    end
end

function value = local_positive_integer(value, name)
    if ~isscalar(value) || ~isreal(value) || ~isfinite(value) || ...
            value < 1 || value ~= round(value)
        error('config_params_2024_v2:InvalidParameter', ...
            '%s must be a positive integer.', name);
    end
end

function local_positive_scalar(value, name)
    if ~isscalar(value) || ~isreal(value) || ~isfinite(value) || value <= 0
        error('config_params_2024_v2:InvalidParameter', ...
            '%s must be a positive real scalar.', name);
    end
end

function local_nonnegative_scalar(value, name)
    if ~isscalar(value) || ~isreal(value) || ~isfinite(value) || value < 0
        error('config_params_2024_v2:InvalidParameter', ...
            '%s must be a nonnegative real scalar.', name);
    end
end