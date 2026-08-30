function params = config_params_2022_v2(overrides)
%CONFIG_PARAMS_2022_V2 Configure the main-physical Srivastava 2022 model.
    if nargin < 1 || isempty(overrides), overrides = struct(); end
    if ~isstruct(overrides) || ~isscalar(overrides)
        error('config_params_2022_v2:InvalidOverrides', ...
            'overrides must be a scalar struct.');
    end
    folder = fileparts(mfilename('fullpath'));
    project_root = fileparts(folder);
    addpath(project_root, '-begin');
    original_folder = pwd;
    cleanup = onCleanup(@() cd(original_folder));
    cd(project_root);
    params = config_params_std(overrides);

    defaults = struct( ...
        'SNR_db_list', -10:5:30, 'trials', 2, 'use_parallel', false, ...
        'verbose_output', true, 'N_pilots_2022', max(1, round(params.N/4)), ...
        'doppler_oversampling_2022', 5, 'data_power_2022', 1, ...
        'pilot_power_2022', 1, 'em_tol_2022', 1e-6, ...
        'em_max_iter_2022', 50, 'variance_floor_2022', 1e-12, ...
        'rgomp_tolerance_2022', 1e-8, 'rgomp_max_groups_2022', 16, ...
        'rgomp_noise_floor_factor_2022', 1.1, ...
        'support_threshold_2022', 1e-3, 'lmmse_tol_2022', 1e-7, ...
        'lmmse_max_iter_2022', 500);
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
        error('config_params_2022_v2:DelayBounds', ...
            'M_tau must cover configured delays and not exceed M.');
    end
    params.delay_grid_2022 = (0:params.M_tau-1).';

    oversampling = local_positive_integer(params.doppler_oversampling_2022, ...
        'doppler_oversampling_2022');
    if isempty(params.doppler_range)
        nu_max = max(abs(configured_dopplers));
    else
        if numel(params.doppler_range) ~= 2 || any(~isfinite(params.doppler_range))
            error('config_params_2022_v2:InvalidDopplerRange', ...
                'doppler_range must be a finite two-element vector.');
        end
        nu_max = max(abs(params.doppler_range));
    end
    required_doppler_count = 2*ceil(nu_max*oversampling)+1;
    if ~isfield(overrides, 'G_nu'), params.G_nu = required_doppler_count; end
    local_positive_integer(params.G_nu, 'G_nu');
    if mod(params.G_nu, 2) ~= 1 || params.G_nu < required_doppler_count
        error('config_params_2022_v2:DopplerBounds', ...
            'G_nu must be odd and cover the configured symmetric Doppler range.');
    end
    center = (params.G_nu-1)/2;
    params.doppler_grid_2022 = ((0:params.G_nu-1)-center).' / oversampling;

    integer_names = {'N_pilots_2022','em_max_iter_2022', ...
        'rgomp_max_groups_2022','lmmse_max_iter_2022'};
    for index = 1:numel(integer_names)
        local_positive_integer(params.(integer_names{index}), integer_names{index});
    end
    if params.N_pilots_2022 > params.N
        error('config_params_2022_v2:PilotCount', ...
            'N_pilots_2022 must not exceed N.');
    end
    local_positive_scalar(params.data_power_2022, 'data_power_2022');
    local_positive_scalar(params.pilot_power_2022, 'pilot_power_2022');
    local_positive_scalar(params.variance_floor_2022, 'variance_floor_2022');
    local_nonnegative_scalar(params.em_tol_2022, 'em_tol_2022');
    local_nonnegative_scalar(params.rgomp_tolerance_2022, 'rgomp_tolerance_2022');
    local_positive_scalar(params.rgomp_noise_floor_factor_2022, ...
        'rgomp_noise_floor_factor_2022');
    local_nonnegative_scalar(params.support_threshold_2022, 'support_threshold_2022');
    local_positive_scalar(params.lmmse_tol_2022, 'lmmse_tol_2022');

    params.padType = 'RCP';
    params.cp_length_2022 = max(params.delay_grid_2022);
    params.total_columns_2022 = params.N_pilots_2022 + params.N;
    params.pilot_locations_2022 = 1:params.N_pilots_2022;
    params.data_locations_2022 = params.N_pilots_2022 + (1:params.N);
    params.pilot_frame_samples_2022 = params.M*params.N_pilots_2022;
    params.data_frame_samples_2022 = params.M*params.N;
    params.pilot_serial_length_2022 = ...
        params.cp_length_2022 + params.pilot_frame_samples_2022;
    params.data_serial_length_2022 = ...
        params.cp_length_2022 + params.data_frame_samples_2022;
    params.doppler_reference_samples_2022 = params.M*params.N;
    params.pilot_overhead_2022 = params.N_pilots_2022 / params.total_columns_2022;
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
            any(~isfinite(dopplers)) || any(delays < 0)
        error('config_params_2022_v2:InvalidTargets', ...
            'Configured delays and Dopplers must be finite and nonempty.');
    end
end

function value = local_positive_integer(value, name)
    if ~isscalar(value) || ~isreal(value) || ~isfinite(value) || ...
            value < 1 || value ~= round(value)
        error('config_params_2022_v2:InvalidParameter', ...
            '%s must be a positive integer.', name);
    end
end

function local_positive_scalar(value, name)
    if ~isscalar(value) || ~isreal(value) || ~isfinite(value) || value <= 0
        error('config_params_2022_v2:InvalidParameter', ...
            '%s must be a positive real scalar.', name);
    end
end

function local_nonnegative_scalar(value, name)
    if ~isscalar(value) || ~isreal(value) || ~isfinite(value) || value < 0
        error('config_params_2022_v2:InvalidParameter', ...
            '%s must be a nonnegative real scalar.', name);
    end
end
