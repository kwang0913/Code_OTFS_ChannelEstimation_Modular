function [Hhat, state] = estimate_rgbl_initial_2024_v2(R, Psi, Rw, params)
%ESTIMATE_RGBL_INITIAL_2024_V2 Run RGBL on the v2 main-channel dictionary.
    [Hhat, state] = estimate_rgbl_initial(R, Psi, Rw, params);
end
