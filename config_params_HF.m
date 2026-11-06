function params = config_params_HF(overrides)
%CONFIG_PARAMS_HF High-frequency (24.25 GHz) simulation parameters.

    params.config_name = 'hf';

    % --- Control ---
    params.show_figures = false;
    params.verbose_output = false;

    params.pulse_shape = 'rectangle';
    params.pilot_power = 12;
    params.pilot_type = 'random';  % 'random' | 'deterministic' | 'zadoff_chu'
    params.pilot_blocks_def = [];  % pilot layout (empty = use defaults in config_pilots_v2)
    params.pilot_block_len = 64;   % delay/Doppler block length (must not exceed N)
    params.n_pilotspertx = 16;     % random pilots per remaining TX antenna (0 = none)
    params.padType = 'RCP';        % 'RCP' | 'CP' | 'ZP' | 'NONE'

    % --- System ---
    params.Nt = 4;
    params.Nr = 16;
    params.fc = 24.25e9;   % Hz
    params.M = 1024;       % delay taps (subcarriers)
    params.N = 64;         % Doppler taps (symbols)
    params.df = 120e3;     % subcarrier spacing (Hz)
    params.c = 299792458;  % exact speed of light in vacuum (m/s)

    % --- Antenna Spacing (in wavelengths) ---
    params.ddlt = 0.5;
    params.ddlr = 0.5;

    params.ModOrder = 4;
    params.modulation = 'QAM';  % 'QAM' | 'PSK'

    % --- Simulation Loop ---
    % params.SNR_db_list = 40;
    params.SNR_db_list = -30:10:40;
    params.trials = 1;
    params.use_parallel = true;    % parfor over trials when trials > 1 (needs Parallel Computing Toolbox)
    params.rng_seed_start = 1;     % trial i uses seed rng_seed_start + i - 1 at every SNR
    params.compute_ber = [];       % [] = detect only when trials == 1; true = every trial

    % --- Target Selection ---
    params.selected_targets = [1, 2, 3, 4];
    params.n_targets = [];         % randomly pick N targets from pool (overrides selected_targets)
    params.targets = [];           % custom targets (empty = use selected_targets from pool)
    params.potential_targets = [];  % custom pool (empty = use defaults)
    params.targets_file = '';       % JSON file for target pool (empty = use built-in default)
    params.target_dopplers = [];    % override Doppler values after selection
    params.doppler_range = [];      % random Doppler with magnitude in [lo, hi], e.g. [3.0, 6.0]
    params.kappa = [];              % fixed fractional Doppler (empty = step by 0.1)

    % --- Tx/Rx Coordinates ---
    params.Coord.Tx.x = 0;
    params.Coord.Tx.y = 0;
    params.Coord.Rx.x = 100;
    params.Coord.Rx.y = 0;
    params.init.R0 = params.Coord.Rx.x - params.Coord.Tx.x;

    % --- SSR Grid: search range (+/- from coarse estimate) ---
    params.ssr.range.a = 8;    % AoA
    params.ssr.range.d = 32;   % AoD
    params.ssr.range.r = 8;    % delay
    params.ssr.range.v = 1;  % Doppler

    % --- SSR Grid: step length ---
    params.ssr.length.a = 1.6;
    params.ssr.length.d = 6.4;
    params.ssr.length.r = 1.6;
    params.ssr.length.v = 0.1;

    % Apply overrides before derived computations
    if nargin > 0 && ~isempty(overrides)
        params = apply_overrides(params, overrides);
    end

    % === Derived ===

    % Pilot config is done per-trial inside generate_TF (RNG-controlled).
    % Pre-compute n_meas/n_pilots_per_tx for sizing only; isolate RNG.
    if ~isfield(params, 'skip_pilot_config') || ~params.skip_pilot_config
        rng_state = rng;
        [~, params.n_meas, params.n_pilots_per_tx] = config_pilots_v2(params);
        rng(rng_state);
    end

    params = resolve_targets(params, 'default_HF.json');
end
