function frame = generate_ap_sip_frame(params)
%GENERATE_AP_SIP_FRAME Create a paper-specific affine-precoded MIMO frame.
%   The N-dimensional unitary split reserves K2 orthogonal columns for
%   deterministic pilots and uses the remaining columns for QAM data.

    local_validate_params(params);

    M = params.M;
    N = params.N;
    Nt = params.Nt;
    K2 = params.K2;
    data_columns = N - K2;
    U = fft(eye(N)) / sqrt(N);

    frame.P = U(:, 1:K2);
    frame.D = U(:, K2 + 1:end);
    frame.source_bits = randi([0, 1], ...
        M * log2(params.ModOrder), data_columns, Nt);
    frame.bits = frame.source_bits;
    frame.Xd = sqrt(local_power(params, 'data_power_2023')) * ...
        local_qam_symbols(frame.source_bits, params.ModOrder);
    frame.Xp = sqrt(local_power(params, 'pilot_power_2023')) * ...
        local_pilot_symbols(M, K2, Nt, params.ModOrder);
    frame.XDD = local_affine_precoding(frame.Xd, frame.D) + ...
        local_affine_precoding(frame.Xp, frame.P);
end

function XDD = local_affine_precoding(symbols, precoder)
    XDD = complex(zeros(size(symbols, 1), size(precoder, 1), size(symbols, 3)));
    for tx = 1:size(symbols, 3)
        XDD(:, :, tx) = symbols(:, :, tx) * precoder';
    end
end

function symbols = local_qam_symbols(bits, mod_order)
    bits_per_symbol = log2(mod_order);
    indices = bit2int(reshape(bits, bits_per_symbol, []), bits_per_symbol);
    symbols = reshape(qammod(indices, mod_order, UnitAveragePower=true), size(bits, 1) / bits_per_symbol, size(bits, 2), size(bits, 3));
end

function pilots = local_pilot_symbols(M, K2, Nt, mod_order)
    pilot_indices = randi([0, mod_order - 1], M, K2, Nt);
    pilots = qammod(pilot_indices, mod_order, UnitAveragePower=true);
end

function local_validate_params(params)
    required = {'M', 'N', 'Nt', 'K2', 'ModOrder'};
    for i = 1:numel(required)
        if ~isfield(params, required{i})
            error('generate_ap_sip_frame:MissingParameter', ...
                'params.%s is required.', required{i});
        end
    end
    integer_fields = {'M', 'N', 'Nt', 'K2'};
    for i = 1:numel(integer_fields)
        value = params.(integer_fields{i});
        if ~isscalar(value) || ~isfinite(value) || value < 1 || value ~= round(value)
            error('generate_ap_sip_frame:InvalidParameter', ...
                'params.%s must be a positive integer.', integer_fields{i});
        end
    end
    if params.K2 >= params.N
        error('generate_ap_sip_frame:InvalidPilotDimension', ...
            'params.K2 must be smaller than params.N.');
    end
    if ~isscalar(params.ModOrder) || params.ModOrder < 4 || ...
            params.ModOrder ~= round(params.ModOrder) || ...
            log2(params.ModOrder) ~= round(log2(params.ModOrder))
        error('generate_ap_sip_frame:InvalidModOrder', ...
            'params.ModOrder must be a power-of-two QAM order of at least four.');
    end
    for name = {'data_power_2023', 'pilot_power_2023'}
        if isfield(params, name{1})
            local_power(params, name{1});
        end
    end
end

function value = local_power(params, name)
    value = 1;
    if isfield(params, name)
        value = params.(name);
    end
    if ~isscalar(value) || ~isreal(value) || ~isfinite(value) || value <= 0
        error('generate_ap_sip_frame:InvalidPower', ...
            'params.%s must be a finite positive real scalar.', name);
    end
end
