function params = config_params_std(overrides)
%CONFIG_PARAMS_STD Standard (4 GHz) simulation parameters.

    params.config_name = 'std';

    % --- Control ---
    params.show_figures = false;
    params.verbose_output = false;

    % --- Pilots ---
    params.pilot_power = 4;
    params.pilot_type = 'random';  % pilot values: 'random' (QPSK) | 'deterministic' (all ones)
    params.pilot_blocks_def = [];  % structured pilot blocks; unless overridden it is default_pilot_blocks(params); [] = random pilots only
    params.n_pilotspertx = 16;     % random pilots per Tx without a block: scalar or one value per Tx (0 = none)
    params.pilot_block_len = 64;   % block length of the default layout (e.g. 32, 64, 128)

    % --- Guard interval ---
    params.padType = 'RCP';        % 'RCP' | 'CP' | 'ZP' | 'RZP' | 'NONE'
    params.padLen = 11;            % guard length (samples), chosen by design; must cover the max total delay of the target pool (delay_range for the automatic grid)

    % --- System ---
    params.Nt = 4;
    params.Nr = 16;
    params.fc = 4e9;       % Hz
    params.M = 512;        % delay taps (subcarriers)
    params.N = 128;        % Doppler taps (symbols)
    params.df = 40e3;      % subcarrier spacing (Hz)
    params.c = 299792458;  % exact speed of light in vacuum (m/s)
    params.range_factor = 1;  % range/velocity resolution factor: 1 bistatic, 2 monostatic round trip
    params.ddlt = 0.5;        % Tx antenna spacing (wavelengths)
    params.ddlr = 0.5;        % Rx antenna spacing (wavelengths)
    params.ModOrder = 4;
    params.modulation = 'QAM';  % 'QAM' | 'PSK'
    % Tx/Rx positions (bistatic geometry)
    params.Coord.Tx.x = 0;
    params.Coord.Tx.y = 0;
    params.Coord.Rx.x = 100;
    params.Coord.Rx.y = 0;

    % --- Simulation loop ---
    params.SNR_db_list = -10:5:30;
    params.trials = 200;
    params.rng_seed_start = 1;     % trial i uses rng(rng_seed_start + i - 1, 'threefry')
    params.n_workers = 6;          % parfor workers when trials > 1; 0 = serial

    % --- Targets: targets (explicit) > n_targets (random) > selected_targets (pool) ---
    % Pool: potential_targets > targets_file > automatic grid (n_targets only) > packaged default.
    params.targets = [];           % explicit DoD/DoA/Doppler/Sign structs
    params.n_targets = [4];        % random targets per trial, drawn from the pool
    params.selected_targets = [1, 6, 11, 16];  % pool indices (targets and n_targets empty)
    params.targets_file = 'targets/pool_std.json';  % JSON pool; empty = automatic grid (n_targets) or packaged default
    params.doppler_range = [3 6];  % random Doppler with magnitude in [lo, hi], e.g. [3.0, 6.0]
    params.kappa = [];             % fixed fractional Doppler (empty = doppler_resolution grid)
    params.aoa_range = [-60 60];   % automatic grid: AoA [lo, hi] (deg, inclusive)
    params.delay_range = [0 11];   % automatic grid: total delay [lo, hi] (samples, inclusive) on the delay_resolution grid
    params.aod_range = [-60 60];   % automatic grid: AoD [lo, hi] (deg, inclusive)
    params.potential_targets = []; % explicit pool (struct array), overrides targets_file
    params.target_dopplers = [];   % override Doppler values after selection
    params.include_los = false;    % static direct path in the channel (known to the SSR, not detected by the coarse stage)

    % --- Resolutions: generation rounds angles/delays/Dopplers to these grids; SSR stops at them by default ---
    params.angle_resolution   = 0.1;  % deg
    params.delay_resolution   = 1;    % samples (integer-delay channel)
    params.doppler_resolution = 0.1;  % bins

    % --- SSR: search half-width (+/- from coarse estimate), step length, min step (binary refinement stops here) ---
    params.ssr.range.a = 8;   % AoA
    params.ssr.range.r = 4;   % delay
    params.ssr.range.v = 2;   % Doppler
    params.ssr.range.d = 64;  % AoD
    params.ssr.length.a = 1.6;
    params.ssr.length.d = 12.8;
    params.ssr.length.r = 0.4;
    params.ssr.length.v = 0.4;
    params.ssr.min_step.a = [];  % [] = angle_resolution
    params.ssr.min_step.r = [];  % [] = delay_resolution
    params.ssr.min_step.v = [];  % [] = doppler_resolution

    % Apply overrides before derived computations
    if nargin > 0 && ~isempty(overrides)
        params = apply_overrides(params, overrides);
    end
    if nargin == 0 || ~isfield(overrides, 'pilot_blocks_def')
        params.pilot_blocks_def = default_pilot_blocks(params);
    end

    validateattributes(params.delay_resolution, {'numeric'}, ...
        {'real','finite','scalar','positive','<=',1}, mfilename, 'delay_resolution');
    ssr_min_steps(params);  % validates params.ssr.min_step

    % === Derived ===
    % n_meas is deterministic from the layout (positions vary, count doesn't)
    rng_state = rng;
    pb = config_pilots(params);
    rng(rng_state);
    params.n_meas = compute_n_meas(pb, params.Nr);

    params.init.R0 = hypot(params.Coord.Rx.x-params.Coord.Tx.x, ...
        params.Coord.Rx.y-params.Coord.Tx.y);

    params = resolve_targets(params, 'pool_std.json');
end
