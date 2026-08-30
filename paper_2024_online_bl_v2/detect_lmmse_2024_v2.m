function [bitsHat, soft, diagnostics] = ...
        detect_lmmse_2024_v2(received, Hhat, frame, params, n0)
%DETECT_LMMSE_2024_V2 Detect data slots with the main CP/Doppler timing.
    local_validate_inputs(received, Hhat, frame, params, n0);
    M = params.M;
    Q = params.M_tau * params.G_nu;
    data_locations = frame.data_locations;
    time_hat = complex(zeros(M, numel(data_locations), params.Nt));
    diagnostics = struct('iterations', zeros(1, numel(data_locations)), ...
        'relative_residual', zeros(1, numel(data_locations)), ...
        'converged', false(1, numel(data_locations)));
    regularization = n0 / params.data_power_2024;

    for data_index = 1:numel(data_locations)
        slot = data_locations(data_index);
        y = reshape(received(:, slot, :), M*params.Nr, 1);
        forward = @(x) local_apply(x, Hhat, params, Q, slot);
        adjoint = @(z) local_adjoint(z, Hhat, params, Q, slot);
        rhs = adjoint(y);
        effective_regularization = regularization;
        if effective_regularization == 0
            effective_regularization = eps(max(1, norm(rhs)));
        end
        dimension = M*params.Nt;
        if dimension <= 2048
            channel = local_channel_block(Hhat, params, Q, slot);
            normal = channel'*channel + ...
                effective_regularization*eye(dimension);
            normal = (normal+normal')/2;
            xhat = normal\rhs;
            diagnostics.iterations(data_index) = 1;
            diagnostics.relative_residual(data_index) = ...
                norm(normal*xhat-rhs)/max(1,norm(rhs));
            diagnostics.converged(data_index) = true;
        else
            normal_map = @(x) adjoint(forward(x)) + ...
                effective_regularization*x;
            [xhat, flag, relative_residual, iterations] = pcg( ...
                normal_map, rhs, local_tolerance(params), ...
                local_max_iterations(params));
            diagnostics.iterations(data_index) = iterations;
            diagnostics.relative_residual(data_index) = relative_residual;
            diagnostics.converged(data_index) = flag == 0;
            if flag ~= 0
                error('detect_lmmse_2024_v2:NotConverged', ...
                    'LMMSE PCG failed with flag %d and residual %.3e.', ...
                    flag, relative_residual);
            end
        end
        time_hat(:, data_index, :) = reshape(xhat, M, 1, params.Nt);
    end

    unitary_dft = fft(eye(numel(data_locations))) / ...
        sqrt(numel(data_locations));
    soft = complex(zeros(size(time_hat)));
    for tx = 1:params.Nt
        soft(:, :, tx) = time_hat(:, :, tx)*unitary_dft;
    end
    hard = qamdemod(soft(:)/sqrt(params.data_power_2024), ...
        params.ModOrder, UnitAveragePower=true);
    bitsHat = reshape(int2bit(hard, log2(params.ModOrder)), ...
        [M*log2(params.ModOrder), numel(data_locations), params.Nt]);
end

function y = local_apply(x, Hhat, params, Q, slot)
    M = params.M;
    y = complex(zeros(M, params.Nr));
    phase = local_phase(params, slot);
    for tx = 1:params.Nt
        input = x((tx-1)*M+(1:M));
        for delay_index = 1:params.M_tau
            rows = (tx-1)*Q + delay_index:params.M_tau:tx*Q;
            gains = phase*Hhat(rows, :);
            y = y + gains.*circshift(input, ...
                params.delay_grid_2024(delay_index));
        end
    end
    y = y(:);
end

function x = local_adjoint(y, Hhat, params, Q, slot)
    M = params.M;
    y = reshape(y, M, params.Nr);
    x = complex(zeros(M*params.Nt, 1));
    phase = local_phase(params, slot);
    for tx = 1:params.Nt
        accumulator = complex(zeros(M, 1));
        for delay_index = 1:params.M_tau
            rows = (tx-1)*Q + delay_index:params.M_tau:tx*Q;
            gains = phase*Hhat(rows, :);
            accumulator = accumulator + circshift( ...
                sum(conj(gains).*y, 2), ...
                -params.delay_grid_2024(delay_index));
        end
        x((tx-1)*M+(1:M)) = accumulator;
    end
end

function channel = local_channel_block(Hhat, params, Q, slot)
    M = params.M;
    channel = complex(zeros(M*params.Nr, M*params.Nt));
    phase = local_phase(params, slot);
    samples = (0:M-1).';
    for tx = 1:params.Nt
        for delay_index = 1:params.M_tau
            rows = (tx-1)*Q + delay_index:params.M_tau:tx*Q;
            gains = phase*Hhat(rows, :);
            source = mod(samples-params.delay_grid_2024(delay_index), M)+1;
            for rx = 1:params.Nr
                output_rows = (rx-1)*M+(1:M);
                input_columns = (tx-1)*M+source.';
                indices = sub2ind(size(channel), output_rows, input_columns);
                channel(indices) = channel(indices) + gains(:, rx).';
            end
        end
    end
end

function phase = local_phase(params, slot)
    block_length = params.M + params.cp_length_2024;
    output_samples = (slot-1)*block_length + ...
        params.cp_length_2024 + (0:params.M-1).';
    phase = exp(1j*2*pi*output_samples * ...
        params.doppler_grid_2024(:).' / ...
        (params.doppler_reference_columns_2024*block_length));
end

function value = local_max_iterations(params)
    value = 200;
    if isfield(params, 'lmmse_max_iter_2024'), value = params.lmmse_max_iter_2024; end
end

function value = local_tolerance(params)
    value = 1e-8;
    if isfield(params, 'lmmse_tol_2024'), value = params.lmmse_tol_2024; end
end

function local_validate_inputs(received, Hhat, frame, params, n0)
    Q = params.M_tau*params.G_nu;
    if ~isnumeric(received) || size(received,1) ~= params.M || ...
            size(received,2) ~= params.total_columns_2024 || ...
            size(received,3) ~= params.Nr || any(~isfinite(received(:)))
        error('detect_lmmse_2024_v2:InvalidReceived', ...
            'received must be finite M-by-(N+Np)-by-Nr.');
    end
    if ~isnumeric(Hhat) || size(Hhat,1) ~= Q*params.Nt || ...
            size(Hhat,2) ~= params.Nr || any(~isfinite(Hhat(:)))
        error('detect_lmmse_2024_v2:InvalidChannel', ...
            'Hhat has invalid dimensions or values.');
    end
    if ~isstruct(frame) || ~isfield(frame,'data_locations') || ...
            numel(frame.data_locations) ~= params.N
        error('detect_lmmse_2024_v2:InvalidFrame', ...
            'frame must contain N data locations.');
    end
    if ~isscalar(n0) || ~isreal(n0) || ~isfinite(n0) || n0 < 0
        error('detect_lmmse_2024_v2:InvalidNoiseVariance', ...
            'n0 must be finite and nonnegative.');
    end
end
