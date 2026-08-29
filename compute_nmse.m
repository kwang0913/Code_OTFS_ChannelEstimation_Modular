function nmse_val = compute_nmse(chanParams, final_ests, params, derivedParams)
%COMPUTE_NMSE NMSE between true and estimated DD channels (via impulse response).
    HDD_true = generate_DD_channel(chanParams, params, derivedParams);
    HDD_est  = generate_DD_channel(final_ests, params, derivedParams);

    numerator = sum(abs(HDD_true - HDD_est).^2, 'all');
    denominator = sum(abs(HDD_true).^2, 'all');
    nmse_val = numerator / denominator;
end
