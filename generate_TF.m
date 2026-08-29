function [XTF, XDD, params, databits, guard_mask] = generate_TF(params)
%GENERATE_TF Build MIMO TF grid with superimposed pilots.
%   [XTF, XDD, params] = generate_TF(params)
%   [XTF, XDD, params, databits, guard_mask] = generate_TF(params)
%
%   databits:   M*log2(ModOrder) x N x Nt — original data bits per Tx
%   guard_mask: M x N x Nt logical — true at DD positions zeroed as guards

    % Pilot placement (consumes RNG for random mode).
    if ~(isfield(params, 'skip_pilot_config') && params.skip_pilot_config)
        [params.pilot_blocks, params.n_meas, params.n_pilots_per_tx] = config_pilots_v2(params);
    end
    pilot_blocks = params.pilot_blocks;

    % Random QAM data grid
    databits = randi([0 1], params.M*log2(params.ModOrder), params.N, params.Nt);
    int_syms = bit2int(reshape(databits, params.M*log2(params.ModOrder), []), log2(params.ModOrder));
    if isfield(params, 'modulation') && strcmpi(params.modulation, 'PSK')
        datagrid = reshape(pskmod(int_syms, params.ModOrder), params.M, params.N, params.Nt);
    else
        datagrid = reshape(qammod(int_syms, params.ModOrder, UnitAveragePower=true), params.M, params.N, params.Nt);
    end

    % Zero out DD bins at pilot/reference positions
    XDD = datagrid;
    guard_mask = false(params.M, params.N, params.Nt);
    for tx_idx = 1:params.Nt
        temp = XDD(:,:,tx_idx);
        temp = temp(:);
        idx = randperm(params.M*params.N, params.n_meas);
        temp(idx) = 0;
        XDD(:,:,tx_idx) = reshape(temp, params.M, params.N);

        mask_lin = false(params.M*params.N, 1);
        mask_lin(idx) = true;
        guard_mask(:,:,tx_idx) = reshape(mask_lin, params.M, params.N);
    end

    % ISFFT: DD→TF (fft along delay, ifft along Doppler, with normalization)
    XTF = otfs_dd_to_tf(XDD, params.M, params.N);

    % Deduplicate shared references: multiple pilot blocks may share one ref symbol
    ref_keys = zeros(length(pilot_blocks), 3);
    for b = 1:length(pilot_blocks)
        ref_keys(b,:) = [pilot_blocks(b).ref_m, pilot_blocks(b).ref_n, pilot_blocks(b).ref_tx];
    end
    [unique_refs, ~, ref_map] = unique(ref_keys, 'rows');

    ref_vals = zeros(size(unique_refs, 1), 1);
    for r = 1:size(unique_refs, 1)
        if strcmpi(params.pilot_type, 'random')
            ref_vals(r) = exp(1j * (pi/4 + randi([0 3]) * pi/2));
        else
            ref_vals(r) = 1;
        end
    end

    % Pass 1: Zero guard zone — clear all pilot/reference TF positions before inserting pilots
    for b = 1:length(pilot_blocks)
        blk = pilot_blocks(b);
        if strcmp(blk.dim, 'random')
            for p = 1:length(blk.m)
                XTF(blk.m(p), blk.n(p), 1:params.Nt) = 0;
            end
        else
            XTF(blk.m, blk.n, 1:params.Nt) = 0;
        end
        XTF(blk.ref_m, blk.ref_n, 1:params.Nt) = 0;
    end

    % Pass 2: Place pilots and references
    for b = 1:length(pilot_blocks)
        blk = pilot_blocks(b);

        if strcmpi(params.pilot_type, 'random')
            pilot_seq = exp(1j * (pi/4 + randi([0 3], 1, blk.n_pilots) * pi/2));
        elseif strcmpi(params.pilot_type, 'zadoff_chu') || strcmpi(params.pilot_type, 'zc')
            pilot_seq = make_zc_sequence(blk.n_pilots, 1);
        else
            pilot_seq = ones(1, blk.n_pilots);
        end

        if strcmp(blk.dim, 'random')
            for p = 1:blk.n_pilots
                XTF(blk.m(p), blk.n(p), blk.tx) = params.pilot_power * pilot_seq(p);
            end
        elseif strcmp(blk.dim, 'delay')
            XTF(blk.m, blk.n, blk.tx) = params.pilot_power * reshape(pilot_seq, length(blk.m), []);
        else
            XTF(blk.m, blk.n, blk.tx) = params.pilot_power * pilot_seq(:).';
        end

        XTF(blk.ref_m, blk.ref_n, blk.ref_tx) = params.pilot_power * ref_vals(ref_map(b));
    end
end
