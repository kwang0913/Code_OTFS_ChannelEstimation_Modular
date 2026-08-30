function channel = reconstruct_dd_channel_2023(H, params)
%RECONSTRUCT_DD_CHANNEL_2023 DD-link energy without MN-by-MN expansion.
%   Every affine-precoded column shares one M-by-M link.  Frobenius energy is
%   therefore N times the sum of link-block energies.  At production scale
%   the energy is accumulated from the exact delay/Doppler atom Gram matrix;
%   small systems additionally retain explicit blocks for numerical oracles.
    Q = params.M_tau * params.G_nu;
    if ~isnumeric(H) || ~isequal(size(H), [Q * params.Nt, params.Nr]) || ...
            any(~isfinite(H(:)))
        error('reconstruct_dd_channel_2023:InvalidCSI', 'H has incompatible dimensions.');
    end
    M = params.M;
    N = params.N;
    samples = (0:M-1).';
    blocks = {};
    retain_blocks = M <= 64;
    if retain_blocks
        blocks = cell(params.Nr, params.Nt);
    end
    block_energy = 0;
    for rx = 1:params.Nr
        for tx = 1:params.Nt
            coefficients = reshape(H((tx-1)*Q+(1:Q), rx), params.M_tau, params.G_nu);
            link_energy = 0;
            for delay_index = 1:params.M_tau
                delay = params.delay_grid_2023(delay_index);
                phase = exp(1j * 2 * pi * (samples-delay) * ...
                    params.doppler_grid_2023(:).' / (M*N));
                doppler_coefficients = coefficients(delay_index, :).';
                link_energy = link_energy + real(doppler_coefficients' * ...
                    (phase' * phase) * doppler_coefficients);
            end
            block_energy = block_energy + link_energy;
            if retain_blocks
                blocks{rx, tx} = local_block(coefficients, params, samples);
            end
        end
    end
    channel = struct('blocks', {blocks}, 'frobenius_squared', N * block_energy, ...
        'per_column_frobenius_squared', block_energy);
end

function block = local_block(coefficients, params, samples)
    M = params.M;
    block = complex(zeros(M));
    for delay_index = 1:params.M_tau
        delay = params.delay_grid_2023(delay_index);
        phase = exp(1j * 2 * pi * (samples-delay) * ...
            params.doppler_grid_2023(:).' / (M*params.N));
        modulation = phase * coefficients(delay_index, :).';
        block = block + modulation .* circshift(eye(M), delay, 1);
    end
end
