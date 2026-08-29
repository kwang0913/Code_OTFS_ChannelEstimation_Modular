function n0_actual = estimate_noise_variance(rx_OTFS, n0, ests)
%ESTIMATE_NOISE_VARIANCE Effective noise variance from thermal + CSI error.
%   n0_actual = estimate_noise_variance(rx_OTFS, n0, ests)

    P_rx = mean(abs(rx_OTFS(:)).^2);
    n0_thermal = P_rx * n0 / (1 + n0);  % P_rx = P_sig + P_noise => P_noise = P_rx * SNR^{-1}/(1+SNR^{-1})
    if isfield(ests, 'residual_nmse')
        n0_csi = P_rx * ests.residual_nmse;
    else
        n0_csi = 0;
    end
    n0_actual = max(n0_thermal, n0_csi);  % worst-case: use whichever dominates (thermal or CSI error)
end
