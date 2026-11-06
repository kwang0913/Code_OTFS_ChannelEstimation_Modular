function params = config_params_std(overrides)
%CONFIG_PARAMS_STD Standard (4 GHz) simulation parameters.

    params.config_name = 'std';

    % --- Control ---
    params.show_figures = true;
    params.verbose_output = true;

    params.pulse_shape = 'rectangle';
    params.pilot_power = 4;
    params.pilot_type = 'random';  % 'random' | 'deterministic' | 'zadoff_chu'
    params.pilot_blocks_def = [];  % pilot layout (empty = use defaults in config_pilots_v2)
    params.pilot_block_len = 64;   % delay/Doppler block length (e.g. 32, 64, 128)
    params.n_pilotspertx = 16;     % random pilots per remaining TX antenna (0 = none)
    params.padType = 'RCP';        % 'RCP' | 'CP' | 'ZP' | 'NONE'

    % --- System ---
    params.Nt = 4;
    params.Nr = 16;
    params.fc = 4e9;       % Hz
    params.M = 512;        % delay taps (subcarriers)
    params.N = 128;        % Doppler taps (symbols)
    params.df = 40e3;      % subcarrier spacing (Hz)
    params.c = 299792458;  % exact speed of light in vacuum (m/s)

    % --- Antenna Spacing (in wavelengths) ---
    params.ddlt = 0.5;
    params.ddlr = 0.5;

    params.ModOrder = 4;
    params.modulation = 'QAM';  % 'QAM' | 'PSK'

    % --- Simulation Loop ---
    % params.SNR_db_list = 15;
    params.SNR_db_list = -10:5:30;
    params.trials = 200;
    params.use_parallel = true;    % parfor over trials when trials > 1 (needs Parallel Computing Toolbox)
    params.rng_seed_start = 1;     % trial i uses seed rng_seed_start + i - 1 at every SNR
    params.compute_ber = [];       % [] = detect only when trials == 1; true = every trial

    % --- Target Selection ---
    params.selected_targets = [1, 6, 11, 16];
    params.n_targets = [4];         % randomly pick N targets from pool (overrides selected_targets)
    params.targets = [];           % custom targets (empty = use selected_targets from pool)
    params.potential_targets = [];  % custom pool (empty = use defaults)
    params.targets_file = '';       % JSON file for target pool (empty = use built-in default)
    params.target_dopplers = [];    % override Doppler values after selection
    params.doppler_range = [3 6];      % random Doppler with magnitude in [lo, hi], e.g. [3.0, 6.0]
    params.kappa = [];              % fixed fractional Doppler (empty = step by 0.1)

    % --- Tx/Rx Coordinates ---
    params.Coord.Tx.x = 0;
    params.Coord.Tx.y = 0;
    params.Coord.Rx.x = 100;
    params.Coord.Rx.y = 0;
    params.init.R0 = params.Coord.Rx.x-params.Coord.Tx.x;

    % --- SSR Grid: search range (half-width in grid bins, centered on coarse estimate) ---
    params.ssr.range.a = 8;   % AoA search half-width
    params.ssr.range.d = 64;  % AoD search half-width
    params.ssr.range.r = 4;   % delay search half-width
    params.ssr.range.v = 2;   % Doppler search half-width

    % --- SSR Grid: step size (fractional grid bins per search step) ---
    params.ssr.length.a = 1.6;
    params.ssr.length.d = 12.8;
    params.ssr.length.r = 0.4;
    params.ssr.length.v = 0.4;

    % Apply overrides before derived computations
    if nargin > 0 && ~isempty(overrides)
        params = apply_overrides(params, overrides);
    end

    % === Derived ===

    % Pilot config is now done per-trial inside generate_TF (RNG-controlled).
    % Pre-compute n_meas/n_pilots_per_tx for sizing only; isolate RNG.
    if ~isfield(params, 'skip_pilot_config') || ~params.skip_pilot_config
        rng_state = rng;
        [~, params.n_meas, params.n_pilots_per_tx] = config_pilots_v2(params);
        rng(rng_state);
    end

    params = resolve_targets(params, 'pool_std.json');
end
