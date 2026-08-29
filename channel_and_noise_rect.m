function [YTF, YDD, rx_OTFS] = channel_and_noise_rect(XTF, chanParams, params, derivedParams, SNR_db)
%CHANNEL_AND_NOISE_RECT Full PHY pipeline: Heisenberg (TF→time) → time-domain channel → AWGN → Wigner (time→TF/DD).

    fsamp = derivedParams.fsamp;

    %% 1. Heisenberg Transform (Mod to time domain)
    tx_OTFS = zeros(derivedParams.numSamps, params.Nt);
    for tx_idx = 1:params.Nt
        temp = ifft(squeeze(XTF(:,:,tx_idx)), [], 1) * sqrt(params.M);
        tx_OTFS(:, tx_idx) = apply_padding(temp, derivedParams.padLen, params.padType);
    end

    %% 2. Time-Domain Channel
    rx_signal_len = derivedParams.numSamps + derivedParams.padLen;
    n_paths = length(chanParams.pathDelays);
    A_OTFS = zeros(rx_signal_len, n_paths);

    for j_idx = 1:n_paths
        temp = struct(...
            'pathGains', chanParams.pathGains(j_idx), ...
            'pathDelays', chanParams.pathDelays(j_idx), ...
            'pathDopplerFreqs', chanParams.pathDopplers(j_idx)/derivedParams.Tsamp, ...
            'pathDopplers', chanParams.pathDopplers(j_idx));
        for tx_idx = 1:params.Nt
            A_OTFS(1:derivedParams.numSamps+temp.pathDelays, j_idx) = ...
                A_OTFS(1:derivedParams.numSamps+temp.pathDelays, j_idx) + ...
                dopplerChannel(tx_OTFS(:,tx_idx), fsamp, temp) .* ...
                exp(-1j*2*pi*(tx_idx-1)*params.ddlt*sind(chanParams.AoDs(j_idx)));  % AoD steering: phase shift per Tx element
        end
    end
    % AoA steering matrix: Nr x P spatial signature times P x T path signals → Nr x T
    rx_OTFS = (exp(1j*2*pi*(0:params.Nr-1)'*params.ddlr*sind(chanParams.AoAs(1:n_paths))) * A_OTFS.').';

    % 'measured' mode: awgn scales noise relative to signal power (needed
    % because pilot_power modifies the actual signal level)
    rx_OTFS = awgn(rx_OTFS, SNR_db, 'measured');

    %% 3. Wigner Transform (Demod to TF/DD domain)
    YTF = zeros(params.M, params.N, params.Nr);
    YDD = zeros(params.M, params.N, params.Nr);
    for rx_idx = 1:params.Nr
        temp = strip_padding(rx_OTFS(1:derivedParams.numSamps, rx_idx), params.M, params.N, ...
            derivedParams.padLen, derivedParams.offset, params.padType);
        YTF(:,:,rx_idx) = fft(temp, [], 1) / sqrt(params.M);
        YDD(:,:,rx_idx) = otfs_tf_to_dd(YTF(:,:,rx_idx), params.M, params.N);
    end
end
