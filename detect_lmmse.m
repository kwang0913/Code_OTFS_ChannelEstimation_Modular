function [ber, XDD_hat, solver_info] = detect_lmmse(rx_OTFS, XTF, params, final_ests, derivedParams, databits, guard_mask, n0)
%DETECT_LMMSE Joint DD-data LMMSE detection via matrix-free augmented LSQR.
%
%   [ber, XDD_hat] = detect_lmmse(rx_OTFS, XTF, params,
%         final_ests, derivedParams, databits, guard_mask, n0)
%
%   Returns XDD_hat (M x N x Nt) soft symbols, with zero DD guard entries.
%   Optional solver_info reports LSQR convergence, residual and iteration count.

    M  = params.M;
    N  = params.N;
    Nt = params.Nt;
    Nr = params.Nr;
    ModOrder = params.ModOrder;
    ddlr = params.ddlr;
    ddlt = params.ddlt;
    padLen  = derivedParams.padLen;
    padType = params.padType;
    numSamps = derivedParams.numSamps;
    bps = log2(ModOrder);

    %% 1. Build G_full with phase-corrected gains
    delay_est   = min(max(round(final_ests.pathDelays), 0), padLen);
    doppler_est = final_ests.pathDopplers;
    gains_est   = final_ests.pathGains(:);
    AoA_est     = final_ests.AoAs;
    AoD_est     = final_ests.AoDs;

    steer_rx = exp(1j*2*pi*ddlr*(0:Nr-1).'*sind(AoA_est(:).'));
    steer_tx = exp(-1j*2*pi*ddlt*(0:Nt-1).'*sind(AoD_est(:).'));

    % Phase correction for getG convention
    phase_corr = compute_phase_correction(doppler_est, delay_est, padLen, padType, M, N);
    gains_corrected = gains_est .* phase_corr;

    G_full = sparse(Nr*numSamps, Nt*numSamps);
    for r = 1:Nr
        for t = 1:Nt
            h_eff = gains_corrected .* steer_rx(r,:).' .* steer_tx(t,:).';
            chanEst_rt.pathDelays   = delay_est(:).';
            chanEst_rt.pathDopplers = doppler_est(:).';
            chanEst_rt.pathGains    = h_eff(:).';
            G_rt = getG(M, N, chanEst_rt, padLen, padType);
            row_range = (r-1)*numSamps + (1:numSamps);
            col_range = (t-1)*numSamps + (1:numSamps);
            G_full(row_range, col_range) = G_rt;
        end
    end

    %% 2. Pilot cancellation: reconstruct pilot-only RX signal, subtract from received
    XTF_pilots = zeros(M, N, Nt);
    reserved_tf = false(M, N);
    pilot_blocks = params.pilot_blocks;
    for b = 1:length(pilot_blocks)
        blk = pilot_blocks(b);
        if strcmp(blk.dim, 'random')
            % Random pilot coordinates are paired, not a Cartesian product.
            for p = 1:length(blk.m)
                XTF_pilots(blk.m(p), blk.n(p), blk.tx) = XTF(blk.m(p), blk.n(p), blk.tx);
                reserved_tf(blk.m(p), blk.n(p)) = true;
            end
        else
            XTF_pilots(blk.m, blk.n, blk.tx) = XTF(blk.m, blk.n, blk.tx);
            reserved_tf(blk.m, blk.n) = true;
        end
        XTF_pilots(blk.ref_m, blk.ref_n, blk.ref_tx) = XTF(blk.ref_m, blk.ref_n, blk.ref_tx);
        reserved_tf(blk.ref_m, blk.ref_n) = true;
    end

    est_chanParams = struct();
    est_chanParams.pathDelays   = delay_est;
    est_chanParams.pathDopplers = doppler_est;
    est_chanParams.pathGains    = gains_est(:).';
    est_chanParams.AoDs         = AoD_est;
    est_chanParams.AoAs         = AoA_est;
    [~, ~, rx_pilots] = channel_and_noise_rect(XTF_pilots, est_chanParams, params, derivedParams, Inf);

    y_rx = reshape(rx_OTFS(1:numSamps, :), [], 1);
    y_pilot = reshape(rx_pilots(1:numSamps, :), [], 1);
    y_data = y_rx - y_pilot;

    %% 3. Noise + CSI error variance estimation
    n0_actual = estimate_noise_variance(rx_OTFS, n0, final_ests);

    %% 4. Joint LMMSE over data symbols: [G*T; sqrt(n0)*I] * d
    % T embeds DD data, applies ISFFT and TF masking, then adds time-domain
    % padding. Its adjoint must add CP contributions back to the frame tail.
    frame = struct('M', M, 'N', N, 'Nt', Nt, 'padLen', padLen, ...
        'padType', padType, 'numSamps', numSamps, ...
        'data_indices', find(~guard_mask), ...
        'reserved_tf', repmat(reserved_tf, 1, 1, Nt));
    n_data = numel(frame.data_indices);
    sqrt_noise = sqrt(n0_actual);
    G_aug = @(v, flag) joint_data_operator(v, flag, G_full, frame, sqrt_noise);
    y_aug = [y_data; zeros(n_data, 1)];

    [data_hat, flag, relres, iterations] = lsqr(G_aug, y_aug, 1e-6, 500);
    solver_info = struct('flag', flag, 'converged', flag == 0, ...
        'relative_residual', relres, 'iterations', iterations);
    if flag ~= 0
        warning('detect_lmmse:LSQRNotConverged', ...
            'Joint LSQR returned flag %d after %d iterations (relative residual %.3g).', ...
            flag, iterations, relres);
    elseif params.verbose_output
        fprintf('Joint LSQR converged after %d iterations (relative residual %.3g).\n', ...
            iterations, relres);
    end

    %% 5. Embed recovered data into the DD grid; guards remain exactly zero
    XDD_hat = zeros(M, N, Nt);
    XDD_hat(frame.data_indices) = data_hat;

    %% 6. BER computation (data positions only)
    total_bit_errors = 0;
    total_data_bits = 0;

    for t = 1:Nt
        x_det = XDD_hat(:,:,t);

        if isfield(params, 'modulation') && strcmpi(params.modulation, 'PSK')
            det_int = pskdemod(x_det, ModOrder);
        else
            det_int = qamdemod(x_det, ModOrder, UnitAveragePower=true);
        end
        det_bits = int2bit(det_int(:), bps);
        det_bits = reshape(det_bits, bps, M*N);

        orig_bits = squeeze(databits(:,:,t));
        orig_bits = reshape(orig_bits, bps, M*N);

        data_mask = ~guard_mask(:,:,t);
        data_mask_lin = data_mask(:).';
        total_bit_errors = total_bit_errors + sum(det_bits(:,data_mask_lin) ~= orig_bits(:,data_mask_lin), 'all');
        total_data_bits = total_data_bits + bps * sum(data_mask_lin);
    end

    ber = total_bit_errors / total_data_bits;
end

function out = joint_data_operator(in, flag, G, frame, sqrt_noise)
    if strcmp(flag, 'notransp')
        out = [G * apply_otfs_data_map(in, 'notransp', frame); sqrt_noise * in];
    else
        n_rx = size(G, 1);
        out = apply_otfs_data_map(G' * in(1:n_rx), 'transp', frame) ...
            + sqrt_noise * in(n_rx+1:end);
    end
end
