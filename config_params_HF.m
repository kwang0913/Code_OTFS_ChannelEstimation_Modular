function params = config_params_HF(overrides)
%CONFIG_PARAMS_HF High-frequency (24.25 GHz) simulation parameters.

    params.config_name = 'hf';

    % --- Control ---
    params.show_figures = false;
    params.verbose_output = false;

    % --- Pilots ---
    params.pilot_power = 12;
    params.pilot_type = 'random';  % pilot values: 'random' (QPSK) | 'deterministic' (all ones)
    params.pilot_blocks_def = [];  % structured pilot blocks; unless overridden it is default_pilot_blocks(params); [] = random pilots only
    params.n_pilotspertx = 16;     % random pilots per Tx without a block: scalar or one value per Tx (0 = none)
    params.pilot_block_len = 64;   % block length of the default layout

    % --- Guard interval ---
    params.padType = 'RCP';        % 'RCP' | 'CP' | 'ZP' | 'RZP' | 'NONE'
    params.padLen = 56;            % guard length (samples), chosen by design; must cover the max total delay of the target pool

    % --- System ---
    params.Nt = 4;
    params.Nr = 16;
    params.fc = 24.25e9;   % Hz
    params.M = 1024;       % delay taps (subcarriers)
    params.N = 64;         % Doppler taps (symbols)
    params.df = 120e3;     % subcarrier spacing (Hz)
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
    params.SNR_db_list = -30:10:40;
    params.trials = 1;
    params.rng_seed_start = 1;     % trial i uses rng(rng_seed_start + i - 1, 'threefry')
    params.n_workers = 6;          % parfor workers when trials > 1; 0 = serial

    % --- Targets ---
    params.targets = [];           % custom targets (empty = use selected_targets from pool)
    params.n_targets = [];         % randomly pick N targets from pool (overrides selected_targets)
    params.selected_targets = [1, 2, 3, 4];
    params.targets_file = '';       % JSON file for target pool (empty = use built-in default)
    params.doppler_range = [];      % random Doppler with magnitude in [lo, hi], e.g. [3.0, 6.0]
    params.kappa = [];              % fixed fractional Doppler (empty = doppler_resolution grid)
    params.potential_targets = [];  % custom pool (empty = use defaults)
    params.target_dopplers = [];    % override Doppler values after selection
    params.include_los = false;    % static direct path in the channel (known to the SSR, not detected by the coarse stage)

    % --- Resolutions: generation rounds angles/delays/Dopplers to these grids; SSR stops at them by default ---
    params.angle_resolution   = 0.1;  % deg
    params.delay_resolution   = 1;    % samples (integer-delay channel)
    params.doppler_resolution = 0.1;  % bins

    % --- SSR: search half-width (+/- from coarse estimate), step length, min step (binary refinement stops here) ---
    params.ssr.range.a = 8;   % AoA
    params.ssr.range.r = 8;   % delay
    params.ssr.range.v = 1;   % Doppler
    params.ssr.range.d = 32;  % AoD
    params.ssr.length.a = 1.6;
    params.ssr.length.d = 6.4;
    params.ssr.length.r = 1.6;
    params.ssr.length.v = 0.1;
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

    params = resolve_targets(params, 'default_HF.json');
end
