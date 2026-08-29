function [pilot_blocks, n_meas_total, n_pilots_per_tx] = config_pilots(params)
%CONFIG_PILOTS Original pilot layout: 4 fixed blocks (2 delay + 2 Doppler), one per Tx antenna.
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
        % --- Hard-coded defaults ---
        b = 0;

        % Delay (range) block
        b = b + 1;
        pb(b).dim    = 'delay';
        pb(b).m      = 256-(64-1):256;
        pb(b).n      = 64;
        pb(b).tx     = 1;
        pb(b).ref_m  = 256-(64-1);
        pb(b).ref_n  = 64;
        pb(b).ref_tx = 1;

        % Delay (range) block
        b = b + 1;
        pb(b).dim    = 'delay';
        pb(b).m      = 257:257+(64-1);
        pb(b).n      = 65;
        pb(b).tx     = 3;
        pb(b).ref_m  = 256-(64-1);
        pb(b).ref_n  = 65;
        pb(b).ref_tx = 3;

        % Doppler (velocity) block
        b = b + 1;
        pb(b).dim    = 'doppler';
        pb(b).m      = 257;
        pb(b).n      = 64-(64-1):64;
        pb(b).tx     = 2;
        pb(b).ref_m  = 257;
        pb(b).ref_n  = 64-(64-1);
        pb(b).ref_tx = 2;

        % Doppler (velocity) block
        b = b + 1;
        pb(b).dim    = 'doppler';
        pb(b).m      = 256;
        pb(b).n      = 65:65+(64-1);
        pb(b).tx     = 4;
        pb(b).ref_m  = 256;
        pb(b).ref_n  = 64-(64-1);
        pb(b).ref_tx = 4;


    end

    % ================================================================
    %  Derived quantities
    % ================================================================
    for i = 1:length(pb)
        if strcmp(pb(i).dim, 'delay')
            pb(i).offsets = pb(i).m - pb(i).ref_m;
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
