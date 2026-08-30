function frame = generate_online_pilot_frame_2024_v2(params)
%GENERATE_ONLINE_PILOT_FRAME_2024_V2 Build the N+Np online-pilot payload.
    local_validate_params(params);
    M = params.M;
    N = params.N;
    Nt = params.Nt;
    Np = params.N_pilots_2024;
    total_columns = N + Np;

    pilot_locations = round(linspace(2, total_columns - 1, Np));
    if numel(unique(pilot_locations)) ~= Np
        error('generate_online_pilot_frame_2024_v2:DuplicatePilotLocations', ...
            'Equally spaced pilot locations must be distinct.');
    end
    data_locations = setdiff(1:total_columns, pilot_locations, 'stable');

    bits_per_symbol = log2(params.ModOrder);
    source_bits = randi([0, 1], M * bits_per_symbol, N, Nt);
    data_dd = sqrt(local_power(params, 'data_power_2024')) * ...
        local_qam_symbols(source_bits, params.ModOrder);
    unitary_dft = fft(eye(N)) / sqrt(N);
    data_time = complex(zeros(M, N, Nt));
    for tx = 1:Nt
        data_time(:, :, tx) = data_dd(:, :, tx) * unitary_dft';
    end

    pilots = sqrt(local_power(params, 'pilot_power_2024')) * ...
        local_pilot_vectors(M, Np, Nt);
    transmitted = complex(zeros(M, total_columns, Nt));
    transmitted(:, data_locations, :) = data_time;
    for pilot = 1:Np
        transmitted(:, pilot_locations(pilot), :) = pilots(:, :, pilot);
    end

    data_mask = false(M, total_columns, Nt);
    data_mask(:, data_locations, :) = true;
    frame = struct( ...
        'source_bits', source_bits, ...
        'bits', source_bits, ...
        'data_dd', data_dd, ...
        'data_time', data_time, ...
        'pilots', pilots, ...
        'pilot_locations', pilot_locations, ...
        'data_locations', data_locations, ...
        'data_mask', data_mask, ...
        'transmitted', transmitted, ...
        'payload', transmitted, ...
        'cp_frame', add_column_cp_2024_v2(transmitted, params));
end

function symbols = local_qam_symbols(bits, mod_order)
    bits_per_symbol = log2(mod_order);
    indices = bit2int(reshape(bits, bits_per_symbol, []), bits_per_symbol);
    symbols = reshape(qammod(indices, mod_order, UnitAveragePower=true), ...
        size(bits, 1) / bits_per_symbol, size(bits, 2), size(bits, 3));
end

function pilots = local_pilot_vectors(M, Np, Nt)
    pilot_indices = randi([0, 3], M, Nt, Np);
    pilots = qammod(pilot_indices, 4, UnitAveragePower=true);
end

function local_validate_params(params)
    required = {'M','N','Nt','N_pilots_2024','ModOrder', ...
        'total_columns_2024','cp_length_2024'};
    for index = 1:numel(required)
        if ~isfield(params, required{index})
            error('generate_online_pilot_frame_2024_v2:MissingParameter', ...
                'params.%s is required.', required{index});
        end
    end
    integer_fields = {'M','N','Nt','N_pilots_2024','total_columns_2024', ...
        'cp_length_2024'};
    for index = 1:numel(integer_fields)
        value = params.(integer_fields{index});
        if ~isscalar(value) || ~isreal(value) || ~isfinite(value) || ...
                value < 0 || value ~= round(value)
            error('generate_online_pilot_frame_2024_v2:InvalidDimension', ...
                'params.%s must be a finite nonnegative integer.', ...
                integer_fields{index});
        end
    end
    if params.M < 1 || params.N < 1 || params.Nt < 1 || ...
            params.N_pilots_2024 < 1 || ...
            params.total_columns_2024 ~= params.N + params.N_pilots_2024
        error('generate_online_pilot_frame_2024_v2:InvalidDimension', ...
            'Frame dimensions are inconsistent.');
    end
    if params.N_pilots_2024 > params.N
        error('generate_online_pilot_frame_2024_v2:TooManyPilots', ...
            'N_pilots_2024 must not exceed N.');
    end
    if ~isscalar(params.ModOrder) || params.ModOrder < 4 || ...
            log2(params.ModOrder) ~= round(log2(params.ModOrder))
        error('generate_online_pilot_frame_2024_v2:InvalidModOrder', ...
            'ModOrder must be a power-of-two QAM order of at least four.');
    end
end

function value = local_power(params, name)
    value = params.(name);
    if ~isscalar(value) || ~isreal(value) || ~isfinite(value) || value <= 0
        error('generate_online_pilot_frame_2024_v2:InvalidPower', ...
            'params.%s must be a finite positive real scalar.', name);
    end
end
