function [ber, XDD_hat] = detect_lmmse(rx_OTFS, XTF, params, final_ests, derivedParams, databits, guard_mask, n0)
%DETECT_LMMSE MIMO LMMSE detection via augmented lsqr (Tikhonov regularization).
%
%   [ber, XDD_hat] = detect_lmmse(rx_OTFS, XTF, params,
%         final_ests, derivedParams, databits, guard_mask, n0)
%
%   Returns XDD_hat (M x N x Nt) hard decisions in DD domain.

    M  = params.M;
    N  = params.N;
    Nt = params.Nt;
    Nr = params.Nr;
    ModOrder = params.ModOrder;
    ddlr = params.ddlr;
    ddlt = params.ddlt;
    padLen  = derivedParams.padLen;
    padType = params.padType;
    offset  = derivedParams.offset;
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
    pilot_blocks = params.pilot_blocks;
    for b = 1:length(pilot_blocks)
        blk = pilot_blocks(b);
        XTF_pilots(blk.m, blk.n, blk.tx) = XTF(blk.m, blk.n, blk.tx);
        XTF_pilots(blk.ref_m, blk.ref_n, blk.ref_tx) = XTF(blk.ref_m, blk.ref_n, blk.ref_tx);
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

    %% 4. LMMSE via lsqr: G_aug = [G; sqrt(n0)*I] implements Tikhonov regularization
    n_tx = Nt * numSamps;
    G_aug = [G_full; sqrt(n0_actual) * speye(n_tx)];
    y_aug = [y_data; zeros(n_tx, 1)];

    if params.verbose_output
        x_hat_vec = lsqr(G_aug, y_aug, 1e-6, 500);
    else
        [x_hat_vec, ~] = lsqr(G_aug, y_aug, 1e-6, 500);
    end
    x_hat = reshape(x_hat_vec, numSamps, Nt);

    %% 5. Demodulate: strip CP -> time-to-freq FFT -> Heisenberg to DD domain
    XDD_hat = zeros(M, N, Nt);
    for t = 1:Nt
        x_stripped = strip_padding(x_hat(:,t), M, N, padLen, offset, padType);
        XTF_hat = fft(x_stripped, [], 1) / sqrt(M);
        XDD_hat(:,:,t) = fft(ifft(XTF_hat, [], 1), [], 2) / sqrt(N/M);
    end

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
