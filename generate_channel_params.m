function [chanParams, derivedParams] = generate_channel_params(params)
%GENERATE_CHANNEL_PARAMS Generate all per-trial channel parameters.
%   Handles target pool sampling (if pool mode), then generates gains,
%   Dopplers, and derives AoD/AoA from bistatic geometry.

    K = params.targets_num;

    % --- Target selection (pool mode vs fixed mode) ---
    if ~isempty(params.n_targets) && isfield(params, 'target_pool')
        % Pool mode: randomly sample targets from pool
        sel = sort(randperm(numel(params.target_pool), params.n_targets));
        tl  = params.target_pool(sel);
        chanParams.selected_targets = sel;
        DoDs       = [tl.DoD];
        DoAs       = [tl.DoA];
        pathDelays = DoDs + DoAs;
        Signs      = [tl.Sign];
    else
        % Fixed mode: use predefined targets from config
        DoDs       = params.init.DoDs;
        DoAs       = params.init.DoAs;
        pathDelays = params.init.pathDelays;
        Signs      = params.init.Signs;
    end

    chanParams.pathDelays = pathDelays;
    chanParams.pathGains  = (randn(1, K) + 1j*randn(1, K)) / sqrt(2);

    % --- Doppler generation ---
    if ~isempty(params.doppler_range)
        lo = params.doppler_range(1);
        hi = params.doppler_range(2);
        if ~isempty(params.kappa)
            candidates = (floor(lo):floor(hi)) + params.kappa;
            candidates = candidates(candidates >= lo & candidates <= hi);
        else
            candidates = lo:0.1:hi;
        end
        mag  = candidates(randi(numel(candidates), 1, K));
        sign = 2*(randi(2, 1, K) - 1) - 1;   % random +1/-1 (approaching vs receding)
        chanParams.pathDopplers = sign .* mag;
    else
        if ~isempty(params.n_targets) && isfield(params, 'target_pool')
            chanParams.pathDopplers = [tl.Doppler];
        else
            chanParams.pathDopplers = params.init.pathDopplers;
        end
    end

    % --- Pad-dependent timing ---
    padLen = max(chanParams.pathDelays);
    [Meff, numSamps, offset] = compute_padding_params(params.M, params.N, padLen, params.padType);

    derivedParams.fsamp = params.M * params.df;

    % --- Compute target coordinates from delay taps and geometry ---
    Coord_tgts.x = zeros(1, K);
    Coord_tgts.y = zeros(1, K);
    for i = 1:K
        [Coord_tgts.x(i), Coord_tgts.y(i)] = compute_target_coord(...
            DoDs(i), DoAs(i), params.Coord, derivedParams.fsamp, params.c);
        Coord_tgts.y(i) = Coord_tgts.y(i) * Signs(i);
    end

    % Derive angles: atand(dy/|dx|) converts Cartesian coords to angle (degrees)
    chanParams.AoDs = round(atand((Coord_tgts.y - params.Coord.Tx.y) ./ abs(Coord_tgts.x - params.Coord.Tx.x)), 1);
    chanParams.AoAs = round(atand((Coord_tgts.y - params.Coord.Rx.y) ./ abs(Coord_tgts.x - params.Coord.Rx.x)), 1);

    derivedParams.padLen = padLen;
    derivedParams.Meff = Meff;
    derivedParams.offset = offset;
    derivedParams.numSamps = numSamps;
    derivedParams.dt = Meff / derivedParams.fsamp;
    derivedParams.Tsamp = params.N * derivedParams.dt;
    derivedParams.T_s = 1 / derivedParams.fsamp;
    derivedParams.lambda = params.c / params.fc;
    derivedParams.range_res = params.c / (params.M * params.df);
    derivedParams.velocity_res = derivedParams.lambda / derivedParams.Tsamp;

    if params.verbose_output
        fprintf('Target Coords (X): %s\n', mat2str(round(Coord_tgts.x, 2)));
        fprintf('Target Coords (Y): %s\n', mat2str(round(Coord_tgts.y, 2)));
        fprintf('Range resolution: %5.2f m\n', derivedParams.range_res);
        fprintf('Velocity resolution: %5.2f m/s\n', derivedParams.velocity_res);
        fprintf('--- Channel & System Parameters ---\n');
        fprintf('True AoA:      %s\n', sprintf('%5.2f ', chanParams.AoAs));
        fprintf('True AoD:      %s\n', sprintf('%5.2f ', chanParams.AoDs));
        fprintf('True Delay:    %s\n', sprintf('%5.2f ', chanParams.pathDelays));
        fprintf('True Doppler: %s\n', sprintf('%5.2f ', chanParams.pathDopplers));
        fprintf('---------------------------------\n');
    end
end
