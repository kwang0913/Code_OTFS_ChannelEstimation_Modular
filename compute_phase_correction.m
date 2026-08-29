function phase_corr = compute_phase_correction(doppler_est, delay_est, padLen, padType, M, N)
%COMPUTE_PHASE_CORRECTION Phase correction matching getG's time-origin convention.
%   Inverts the exp(-j*2*pi/MN * v*(l+padLen)) phase from getG so that gains are real-valued.
%   phase_corr = compute_phase_correction(doppler_est, delay_est, padLen, padType, M, N)

    [Meff, ~, ~] = compute_padding_params(M, N, padLen, padType);
    MN_G = Meff * N;

    if strcmpi(padType, 'CP') || strcmpi(padType, 'RCP')
        phase_corr = exp(1j * 2*pi / MN_G * doppler_est(:) .* (delay_est(:) + padLen));
    else
        phase_corr = exp(1j * 2*pi / MN_G * doppler_est(:) .* delay_est(:));
    end
end
