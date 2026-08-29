function Phi_est = build_response_matrix(final_ests, params)
%BUILD_RESPONSE_MATRIX Channel dictionary using full xi_j[n]*H^j[n,m] with fractional Doppler.

    pilot_blocks = params.pilot_blocks;
    Nr = params.Nr;
    M = params.M;
    N = params.N;
    K = length(final_ests.AoAs);

    n_meas = compute_n_meas(pilot_blocks, Nr);

    Phi_est = zeros(n_meas, K);

    % AoA steering: Nr x K
    a_steer = exp(1j*2*pi*params.ddlr*(0:Nr-1).'*sind(final_ests.AoAs(:)'));

    nu_arr    = final_ests.pathDopplers(:)';
    k_arr     = floor(nu_arr);       % integer Doppler index
    kappa_arr = nu_arr - k_arr;      % fractional Doppler remainder in [0,1)
    l_arr     = final_ests.pathDelays(:)';

    idx = 0;
    for b = 1:length(pilot_blocks)
        blk = pilot_blocks(b);

        d_ref   = exp(-1j*2*pi*params.ddlt*(blk.ref_tx - 1)*sind(final_ests.AoDs(:)'));  % AoD steering for reference Tx antenna
        d_pilot = exp(-1j*2*pi*params.ddlt*(blk.tx - 1)*sind(final_ests.AoDs(:)'));  % AoD steering for pilot Tx antenna

        % Reference row
        n0 = blk.ref_n - 1;
        m0 = blk.ref_m - 1;
        xi_H_ref = compute_xi_H(nu_arr, k_arr, kappa_arr, l_arr, n0, m0, N, M);
        Phi_est(idx+(1:Nr), :) = a_steer .* d_ref .* xi_H_ref;
        idx = idx + Nr;

        % Pilot rows
        for p = 1:blk.n_pilots
            if strcmp(blk.dim, 'random')
                m0 = blk.m(p) - 1;
                n0 = blk.n(p) - 1;
            elseif strcmp(blk.dim, 'delay')
                n_m = length(blk.m);
                m_idx = mod(p-1, n_m) + 1;      % column-major: m varies fastest
                n_idx = floor((p-1) / n_m) + 1;  % n varies slowest (blocks of m)
                m0 = blk.m(m_idx) - 1;
                n0 = blk.n(n_idx) - 1;
            else
                n0 = blk.n(p) - 1;
                m0 = blk.m - 1;
            end
            xi_H_p = compute_xi_H(nu_arr, k_arr, kappa_arr, l_arr, n0, m0, N, M);
            Phi_est(idx+(1:Nr), :) = a_steer .* d_pilot .* xi_H_p;
            idx = idx + Nr;
        end
    end
end
