function chanParams = generate_channel_params(params)
%GENERATE_CHANNEL_PARAMS Geometric reflections followed by a static direct path.
%   targets_num counts reflections only. Random target selection samples
%   without replacement. No extra delay is added to any target source.
%   Total delays are quantized using params.delay_resolution before geometry;
%   AoA/AoD and Dopplers are rounded to params.angle_resolution /
%   params.doppler_resolution (an exact params.kappa is kept).
%   An optional static LoS (params.include_los, default false) is LAST, with an unknown
%   gain and no power normalization across paths (the SSR keeps it as a known path; the
%   coarse stage does not detect it). The guard length is fixed by params.padLen; the largest delay must fit in it.
%
%   RNG draw order: targets (selection) -> gains (reflections and LoS) ->
%   Dopplers (magnitude, sign); pilot positions, data, guards, pilot values
%   and noise follow.
    include_los = params.include_los;
    K = params.targets_num;
    use_pool = isfield(params,'target_pool') && ...
        (~isfield(params,'targets') || isempty(params.targets));
    sample_pool = use_pool && ~isempty(params.n_targets);
    if sample_pool
        sel = sort(randperm(numel(params.target_pool), K));
        tl = params.target_pool(sel);
        chanParams.selected_targets = sel;
        DoDs = [tl.DoD];
        DoAs = [tl.DoA];
        Signs = [tl.Sign];
    else
        DoDs = params.init.DoDs;
        DoAs = params.init.DoAs;
        Signs = params.init.Signs;
        if use_pool, chanParams.selected_targets = params.selected_targets; end
    end
    DoDs = DoDs(:).'; DoAs = DoAs(:).'; Signs = Signs(:).';
    total_delay = DoDs + DoAs;
    chanParams.pathDelays = round(total_delay/params.delay_resolution)*params.delay_resolution;
    rounding_error = chanParams.pathDelays-total_delay;
    % Preserve existing grid geometry: decimal addition alone can differ
    % from its rounded value by a few ulps. Only adjust genuinely off-grid totals.
    rounding_error(abs(rounding_error) <= 4*eps(max(1,abs(total_delay)))) = 0;
    delay_adjustment = rounding_error/2;
    DoDs = DoDs + delay_adjustment;
    DoAs = DoAs + delay_adjustment;
    % gains = exp(1j*2*pi*rand(1, K+include_los));                                  % unit modulus, random phase
    gains = sqrt(1/2) * (randn(1, K+include_los) + 1j*randn(1, K+include_los));  % complex Rayleigh

    if ~isempty(params.doppler_range)
        lo = params.doppler_range(1);
        hi = params.doppler_range(2);
        if ~isempty(params.kappa)
            candidates = (floor(lo):floor(hi)) + params.kappa;
            candidates = candidates(candidates >= lo & candidates <= hi);
        else
            candidates = lo:params.doppler_resolution:hi;
        end
        mag = candidates(randi(numel(candidates),1,K));
        signs = 2*(randi(2,1,K)-1)-1;
        chanParams.pathDopplers = signs.*mag;
    elseif sample_pool && isempty(params.target_dopplers)
        chanParams.pathDopplers = [tl.Doppler];
    else
        chanParams.pathDopplers = params.init.pathDopplers;
    end
    chanParams.pathDopplers = chanParams.pathDopplers(:).';
    if isempty(params.kappa)
        chanParams.pathDopplers = snap(chanParams.pathDopplers, params.doppler_resolution);
    end

    fsamp = params.M*params.df;
    tx = [params.Coord.Tx.x; params.Coord.Tx.y];
    rx = [params.Coord.Rx.x; params.Coord.Rx.y];
    baseline = norm(rx-tx);
    baseline_samples = baseline*fsamp/params.c;
    % A target that is feasible before quantization can become infeasible
    % when its total delay rounds below the direct-path delay. Say so instead
    % of reporting a generic triangle error (nothing is moved silently).
    if K > 0 && baseline > 0 && any(total_delay >= baseline_samples & ...
            chanParams.pathDelays < baseline_samples)
        error('generate_channel_params:QuantizedDelayBelowBaseline', ...
            ['Total delay %s rounds (delay_resolution = %g) below the direct-path ' ...
             'delay %.4f samples; use a finer delay_resolution or move the target.'], ...
            mat2str(total_delay(total_delay >= baseline_samples & ...
            chanParams.pathDelays < baseline_samples)), params.delay_resolution, baseline_samples);
    end
    if K > 0 && (baseline == 0 || any(~isfinite(DoDs) | ~isfinite(DoAs) | ...
            DoDs < 0 | DoAs < 0 | DoDs+DoAs < baseline_samples | ...
            abs(DoDs-DoAs) > baseline_samples))
        error('generate_channel_params:InvalidGeometry', ...
            'Reflected range legs must form a triangle with distinct Tx/Rx endpoints.');
    end
    % Reuse the circle-intersection geometry in a baseline-aligned frame,
    % then rotate/translate it; Sign selects a side of the Tx/Rx baseline.
    local.Tx.x = 0; local.Tx.y = 0;
    local.Rx.x = baseline; local.Rx.y = 0;
    Coord_tgts = zeros(2,K);
    if K > 0
        baseline_dir = (rx-tx)/baseline;
        basis = [baseline_dir, [-baseline_dir(2); baseline_dir(1)]];
        for i = 1:K
            [x,y] = compute_target_coord(DoDs(i),DoAs(i),local,fsamp,params.c);
            Coord_tgts(:,i) = tx + basis*[x; Signs(i)*y];
        end
    end
    chanParams.AoDs = snap(atan2d(Coord_tgts(2,:)-tx(2),abs(Coord_tgts(1,:)-tx(1))), params.angle_resolution);
    chanParams.AoAs = snap(atan2d(Coord_tgts(2,:)-rx(2),abs(Coord_tgts(1,:)-rx(1))), params.angle_resolution);

    if include_los
        % The receiver's known LoS uses this same scenario-grid value.
        los = known_los_path(params);
        chanParams.pathDelays = [chanParams.pathDelays, los.pathDelays];
        chanParams.pathDopplers = [chanParams.pathDopplers, los.pathDopplers];
        chanParams.AoDs = [chanParams.AoDs, los.AoDs];
        chanParams.AoAs = [chanParams.AoAs, los.AoAs];
        chanParams.is_los = [false(1,K), true];
    else
        chanParams.is_los = false(1,K);
    end
    chanParams.pathGains = gains;
    if ceil(max(chanParams.pathDelays)) > params.padLen
        error('generate_channel_params:padLen', ...
            'ceil(max delay) = %d exceeds params.padLen = %d.', ceil(max(chanParams.pathDelays)), params.padLen);
    end

    if params.verbose_output
        derivedParams = compute_derived_params(params);
        fprintf('Target Coords (X): %s\n',mat2str(round(Coord_tgts(1,:),2)));
        fprintf('Target Coords (Y): %s\n',mat2str(round(Coord_tgts(2,:),2)));
        fprintf('Range resolution: %5.2f m\n',derivedParams.range_res);
        fprintf('Velocity resolution: %5.2f m/s\n',derivedParams.velocity_res);
        if include_los
            fprintf('--- Channel & System Parameters (LoS last) ---\n');
        else
            fprintf('--- Channel & System Parameters (reflections only) ---\n');
        end
        fprintf('True AoA:     %s\n',sprintf('%5.2f ',chanParams.AoAs));
        fprintf('True AoD:     %s\n',sprintf('%5.2f ',chanParams.AoDs));
        fprintf('True Delay:   %s\n',sprintf('%5.2f ',chanParams.pathDelays));
        fprintf('True Doppler: %s\n',sprintf('%5.2f ',chanParams.pathDopplers));
        fprintf('--------------------------------------------\n');
    end
end

function x = snap(x, res)
% Round to the resolution grid (12-digit rounding removes multiplication noise).
    x = round(round(x / res) * res, 12);
end
