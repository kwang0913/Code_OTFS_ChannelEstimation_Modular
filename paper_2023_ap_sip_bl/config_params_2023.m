function params = config_params_2023(overrides)
%CONFIG_PARAMS_2023 Configure the 2023 AP-SIP sparse MIMO channel model.

    if nargin < 1 || isempty(overrides)
        overrides = struct();
    end
    if ~isstruct(overrides)
        error('config_params_2023:InvalidOverrides', ...
            'overrides must be a struct.');
    end

    paper_dir = fileparts(mfilename('fullpath'));
    project_root = fileparts(paper_dir);
    addpath(project_root, '-begin');
    original_folder = pwd;
    folder_cleanup = onCleanup(@() cd(original_folder));
    cd(project_root);

    % The standard configuration owns physical and target selection settings.
    % Apply paper-specific defaults afterwards so root derived computations do
    % not overwrite paper fields supplied by the caller.
    % The AP-SIP frame is independent of the root delay-Doppler pilot layout,
    % so the root layout is only sized with random pilots unless overridden.
    base_overrides = overrides;
    if ~isfield(overrides, 'pilot_blocks_def')
        base_overrides.pilot_blocks_def = [];
    end
    params = config_params_std(base_overrides);
    params.config_name = 'ap_sip_bl_2023';
    clear folder_cleanup;
    paper_defaults = struct( ...
        'SNR_db_list', -10:5:30, 'trials', 30, 'K2', 8, ...
        'em_tol_2023', 1e-6, 'em_max_iter_2023', 50, ...
        'da_max_iter_2023', 20, 'variance_floor_2023', 1e-12, ...
        'paper_gain_mode', 'fading', 'padType', 'CP', ...
        'use_parallel', false, 'compute_ber', [], ...
        'doppler_oversampling_2023', 4, ...
        'data_power_2023', 0.5, 'pilot_power_2023', 0.5, ...
        'snr_reference_2023', 'data');
    default_fields = fieldnames(paper_defaults);
    for i = 1:numel(default_fields)
        field = default_fields{i};
        if ~isfield(overrides, field)
            params.(field) = paper_defaults.(field);
        end
    end
    if isfield(overrides, 'pilot_power') && ...
            ~isfield(overrides, 'pilot_power_2023') && ~isfield(overrides, 'data_power_2023')
        inherited_pilot_power = overrides.pilot_power;
        if ~isscalar(inherited_pilot_power) || ~isreal(inherited_pilot_power) || ...
                ~isfinite(inherited_pilot_power) || inherited_pilot_power <= 0 || ...
                inherited_pilot_power >= 1
            error('config_params_2023:InvalidPower', ...
                'root pilot_power must be a real scalar strictly between zero and one.');
        end
        params.pilot_power_2023 = inherited_pilot_power;
        params.data_power_2023 = 1 - inherited_pilot_power;
    end

    local_validate_power_contract(params, '2023');

    if ~isscalar(params.doppler_oversampling_2023) || ...
            ~isfinite(params.doppler_oversampling_2023) || ...
            params.doppler_oversampling_2023 < 1 || ...
            params.doppler_oversampling_2023 ~= round(params.doppler_oversampling_2023)
        error('config_params_2023:InvalidDopplerOversampling', ...
            'doppler_oversampling_2023 must be a positive integer.');
    end

    [configured_delays, configured_dods, configured_doas] = local_configured_targets(params);
    required_M_tau = ceil(max(configured_delays)) + 1;
    if ~isfield(overrides, 'M_tau')
        params.M_tau = required_M_tau;
    end
    if ~isscalar(params.M_tau) || ~isfinite(params.M_tau) || ...
            params.M_tau < 1 || params.M_tau ~= round(params.M_tau) || ...
            params.M_tau < required_M_tau
        error('config_params_2023:DelayBounds', ...
            'M_tau must cover every configured target delay plus one guard tap.');
    end
    if strcmpi(params.padType, 'CP') && params.M_tau - 1 > params.M
        error('config_params_2023:CPTooLong', ...
            'Per-symbol CP length M_tau-1 must not exceed M.');
    end
    params.delay_grid_2023 = (0:params.M_tau-1).';

    configured_dopplers = local_configured_dopplers(params);
    if isempty(params.doppler_range)
        nu_max = max(abs(configured_dopplers));
    else
        if numel(params.doppler_range) ~= 2 || ...
                any(~isfinite(params.doppler_range)) || ...
                any(params.doppler_range < 0) || ...
                params.doppler_range(1) > params.doppler_range(2)
            error('config_params_2023:InvalidDopplerRange', ...
                'doppler_range must be an ascending, nonnegative two-element vector.');
        end
        nu_max = max(abs(params.doppler_range));
    end
    oversampling = params.doppler_oversampling_2023;
    required_G_nu = 2 * ceil(nu_max * oversampling) + 1;
    if ~isfield(overrides, 'G_nu')
        params.G_nu = required_G_nu;
    end
    if ~isscalar(params.G_nu) || ~isfinite(params.G_nu) || ...
            params.G_nu < 1 || params.G_nu ~= round(params.G_nu) || ...
            mod(params.G_nu, 2) ~= 1 || params.G_nu < required_G_nu
        error('config_params_2023:DopplerBounds', ...
            'G_nu must be odd and cover the symmetric configured Doppler range.');
    end
    center = (params.G_nu - 1) / 2;
    params.doppler_grid_2023 = ((0:params.G_nu-1) - center).' / oversampling;
    if min(params.doppler_grid_2023) > -nu_max || ...
            max(params.doppler_grid_2023) < nu_max
        error('config_params_2023:DopplerBounds', ...
            'G_nu does not cover the symmetric configured Doppler range.');
    end
end

function [delays, dods, doas] = local_configured_targets(params)
    if ~isempty(params.n_targets) && isfield(params, 'target_pool')
        source = params.target_pool;
        dods = [source.DoD];
        doas = [source.DoA];
    else
        dods = params.init.DoDs;
        doas = params.init.DoAs;
    end
    delays = dods + doas;
    if isempty(delays) || any(~isfinite(delays)) || any(delays < 0)
        error('config_params_2023:InvalidTargetDelays', ...
            'Configured target delays must be finite and nonnegative.');
    end
end


function dopplers = local_configured_dopplers(params)
    if ~isempty(params.n_targets) && isfield(params, 'target_pool')
        source = params.target_pool;
        dopplers = [source.Doppler];
    else
        dopplers = params.init.pathDopplers;
    end
    if isempty(dopplers) || any(~isfinite(dopplers))
        error('config_params_2023:InvalidTargetDopplers', ...
            'Configured target Dopplers must be finite.');
    end
end

function local_validate_power_contract(params, paper)
    fields = {['data_power_' paper], ['pilot_power_' paper]};
    for i = 1:numel(fields)
        value = params.(fields{i});
        if ~isscalar(value) || ~isreal(value) || ~isfinite(value) || value <= 0
            error('config_params_2023:InvalidPower', '%s must be a positive real scalar.', fields{i});
        end
    end
    if abs(params.data_power_2023 + params.pilot_power_2023 - 1) > 100 * eps
        error('config_params_2023:PowerNormalization', ...
            'data_power_2023 + pilot_power_2023 must equal one.');
    end
    reference = params.(['snr_reference_' paper]);
    if ~(ischar(reference) || (isstring(reference) && isscalar(reference))) || ...
            ~strcmpi(char(reference), 'data')
        error('config_params_2023:InvalidSNRReference', ...
            'snr_reference_2023 must be ''data'' (n0 = data_power_2023 / SNR).');
    end
end
