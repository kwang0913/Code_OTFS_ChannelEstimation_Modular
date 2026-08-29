function [pilot_blocks, n_meas_total, n_pilots_per_tx] = config_pilots_v2(params)
%CONFIG_PILOTS_V2 Pilot layout v2: 2 structured blocks (delay+Doppler) + random pilots on remaining Tx antennas.
%
%   Fields per block:
%     dim    - 'delay' or 'doppler' (which dimension the pilots span)
%     m, n   - subcarrier/time-slot indices (array in the spanning dim)
%     tx     - Tx antenna for pilots (1-based)
%     ref_m, ref_n, ref_tx - reference position and antenna
%
%   Derived automatically: offsets, n_pilots
%
%   If params.pilot_blocks_def is set (e.g. from JSON), it is used instead
%   of the hard-coded defaults below.  The m and n fields accept MATLAB
%   range strings such as "2:128" as well as numeric scalars or arrays.

    % ================================================================
    %  PILOT BLOCK DEFINITIONS
    % ================================================================
    if isfield(params, 'pilot_blocks_def') && ~isempty(params.pilot_blocks_def)
        % --- From JSON / override ---
        defs = params.pilot_blocks_def;
        if iscell(defs)
            % jsondecode returns a cell array when struct fields differ
            for i = numel(defs):-1:1
                pb(i) = parse_pilot_block_def(defs{i});
            end
        else
            % struct array (all blocks share the same field set)
            for i = numel(defs):-1:1
                pb(i) = parse_pilot_block_def(defs(i));
            end
        end
    else
        % --- Defaults from params ---
        n_pilotspertx = params.n_pilotspertx;
        blen = params.pilot_block_len;
        b = 0;

        % Delay (range) block, centered in the delay dimension.
        delay_center = floor(params.M / 2);
        if blen > delay_center
            error('config_pilots_v2:blockTooLong', ...
                'pilot_block_len must not exceed floor(M/2).');
        end
        b = b + 1;
        pb(b).dim    = 'delay';
        pb(b).m      = delay_center-(blen-1):delay_center;
        pb(b).n      = floor(params.N / 2);
        pb(b).tx     = 1;
        pb(b).ref_m  = delay_center-(blen-1);
        pb(b).ref_n  = floor(params.N / 2);
        pb(b).ref_tx = 1;

        % Doppler (velocity) block, ending at the final Doppler bin.
        if params.Nt < 2
            error('config_pilots_v2:insufficientTransmitters', ...
                'The structured v2 pilot layout requires Nt >= 2.');
        end
        if blen > params.N
            error('config_pilots_v2:blockTooLong', ...
                'pilot_block_len must not exceed N.');
        end
        b = b + 1;
        pb(b).dim    = 'doppler';
        pb(b).m      = delay_center + 1;
        pb(b).n      = params.N-(blen-1):params.N;
        pb(b).tx     = 2;
        pb(b).ref_m  = delay_center + 1;
        pb(b).ref_n  = params.N-(blen-1);
        pb(b).ref_tx = 2;

        % Random pilots on remaining antennas (tx = 3:Nt)
        M = params.M;
        N = params.N;

        % Collect all [m,n] positions used by blocks above
        used = false(M, N);
        for i = 1:b
            blk_i = pb(i);
            if strcmp(blk_i.dim, 'delay')
                for mi = 1:length(blk_i.m)
                    used(blk_i.m(mi), blk_i.n) = true;
                end
            else  % doppler
                for ni = 1:length(blk_i.n)
                    used(blk_i.m, blk_i.n(ni)) = true;
                end
            end
            used(blk_i.ref_m, blk_i.ref_n) = true;
        end

        [avail_m, avail_n] = find(~used);  % all DD grid positions not occupied by structured blocks
        perm = randperm(length(avail_m));
        avail_m = avail_m(perm);
        avail_n = avail_n(perm);

        if n_pilotspertx > 0 && params.Nt > 2  % assign random pilots to Tx antennas 3..Nt
            offset_pos = 0;
            for tx = 3:params.Nt
                sel = offset_pos + (1:n_pilotspertx);
                b = b + 1;
                pb(b).dim    = 'random';
                pb(b).m      = avail_m(sel)';
                pb(b).n      = avail_n(sel)';
                pb(b).tx     = tx;
                pb(b).ref_m  = avail_m(sel(1));
                pb(b).ref_n  = avail_n(sel(1));
                pb(b).ref_tx = tx;
                offset_pos = offset_pos + n_pilotspertx;
            end
        end

    end

    % ================================================================
    %  Derived quantities
    % ================================================================
    for i = 1:length(pb)
        if strcmp(pb(i).dim, 'delay')
            pb(i).offsets = pb(i).m - pb(i).ref_m;
        elseif strcmp(pb(i).dim, 'random')
            pb(i).offsets = 0:(length(pb(i).m) - 1);
        else
            pb(i).offsets = pb(i).n - pb(i).ref_n;
        end
        pb(i).n_pilots = length(pb(i).offsets);
    end

    % Deduplicated reference count
    ref_keys = zeros(length(pb), 3);
    for i = 1:length(pb)
        ref_keys(i,:) = [pb(i).ref_m, pb(i).ref_n, pb(i).ref_tx];
    end
    unique_refs = unique(ref_keys, 'rows');
    n_unique_refs = size(unique_refs, 1);
    n_meas_total = sum([pb.n_pilots]) + n_unique_refs;

    % Pilots per Tx antenna
    n_pilots_per_tx = zeros(1, params.Nt);
    for i = 1:length(pb)
        n_pilots_per_tx(pb(i).tx) = n_pilots_per_tx(pb(i).tx) + pb(i).n_pilots;
    end
    for i = 1:n_unique_refs
        n_pilots_per_tx(unique_refs(i, 3)) = n_pilots_per_tx(unique_refs(i, 3)) + 1;
    end

    pilot_blocks = pb;
end

% -----------------------------------------------------------------
function blk = parse_pilot_block_def(def)
%PARSE_PILOT_BLOCK_DEF Convert a JSON-decoded struct to a pilot block.
    blk.dim    = char(def.dim);
    blk.m      = parse_range(def.m);
    blk.n      = parse_range(def.n);
    blk.tx     = def.tx;
    blk.ref_m  = parse_range(def.ref_m);
    blk.ref_n  = parse_range(def.ref_n);
    blk.ref_tx = def.ref_tx;
end
