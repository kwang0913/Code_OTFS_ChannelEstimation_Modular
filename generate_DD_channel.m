function HDD_true = generate_DD_channel(chanParams, params, derivedParams)
%GENERATE_DD_CHANNEL True DD channel via per-Tx impulse through noiseless PHY.
%   Sends a single-tap impulse per Tx antenna and collects the Rx DD response.

    % Use the same transmitted waveform for truth and estimates: fractional
    % delays and supplied padding/timing metadata are never recomputed.

    HDD_true = zeros(params.M, params.N, params.Nr, params.Nt);
    for tx_idx = 1:params.Nt
        XDD_pilot = zeros(params.M, params.N, params.Nt);
        XDD_pilot(1, 1, tx_idx) = 1;          % unit impulse at (delay=0, Doppler=0)

        XTF_pilot = otfs_dd_to_tf(XDD_pilot, params.M, params.N);  % ISFFT: DD→TF

        [~, HDD_true(:,:,:,tx_idx), ~] = channel_and_noise_rect(...
            XTF_pilot, chanParams, params, derivedParams, Inf);  % SNR=Inf → noiseless
    end
end
