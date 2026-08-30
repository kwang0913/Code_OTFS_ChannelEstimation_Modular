function nmse = compute_main_nmse_2024_v2( ...
        chanParams, Hhat, params, data_locations)
%COMPUTE_MAIN_NMSE_2024_V2 Main-channel operator NMSE on augmented data slots.
    local_validate_inputs(chanParams, Hhat, params, data_locations);
    block_length = params.M + params.cp_length_2024;
    denominator = params.doppler_reference_columns_2024 * block_length;
    Q = params.M_tau * params.G_nu;
    true_energy = 0;
    error_energy = 0;
    physical_delays = chanParams.pathDelays(:);
    all_delays = unique([params.delay_grid_2024(:); physical_delays]);

    for slot = data_locations(:).'
        output_samples = (slot-1)*block_length + ...
            params.cp_length_2024 + (0:params.M-1).';
        grid_phase = exp(1j*2*pi * output_samples * ...
            params.doppler_grid_2024(:).' / denominator);
        for tx = 1:params.Nt
            for rx = 1:params.Nr
                for delay = all_delays.'
                    true_gain = complex(zeros(params.M, 1));
                    active = find(physical_delays == delay);
                    for index = active.'
                        coefficient = chanParams.pathGains(index) * ...
                            exp(1j*2*pi*(rx-1)*params.ddlr * ...
                                sind(chanParams.AoAs(index))) * ...
                            exp(-1j*2*pi*(tx-1)*params.ddlt * ...
                                sind(chanParams.AoDs(index)));
                        phase = exp(1j*2*pi * output_samples * ...
                            chanParams.pathDopplers(index) / denominator);
                        true_gain = true_gain + coefficient * phase;
                    end

                    estimated_gain = complex(zeros(params.M, 1));
                    delay_index = find(params.delay_grid_2024 == delay, 1);
                    if ~isempty(delay_index)
                        rows = (tx-1)*Q + delay_index:params.M_tau:tx*Q;
                        estimated_gain = grid_phase * Hhat(rows, rx);
                    end
                    true_energy = true_energy + sum(abs(true_gain).^2);
                    error_energy = error_energy + ...
                        sum(abs(estimated_gain-true_gain).^2);
                end
            end
        end
    end
    nmse = real(error_energy) / max(real(true_energy), eps);
    if ~isfinite(nmse) || nmse < 0
        error('compute_main_nmse_2024_v2:NonFiniteMetric', ...
            'Main physical NMSE must be finite and nonnegative.');
    end
end

function local_validate_inputs(chanParams, Hhat, params, data_locations)
    required_params = {'M','N','Nt','Nr','M_tau','G_nu','ddlt','ddlr', ...
        'delay_grid_2024','doppler_grid_2024','cp_length_2024', ...
        'total_columns_2024','doppler_reference_columns_2024'};
    for index = 1:numel(required_params)
        if ~isfield(params, required_params{index})
            error('compute_main_nmse_2024_v2:MissingParameter', ...
                'params.%s is required.', required_params{index});
        end
    end
    Q = params.M_tau * params.G_nu;
    if ~isnumeric(Hhat) || size(Hhat, 1) ~= Q*params.Nt || ...
            size(Hhat, 2) ~= params.Nr || any(~isfinite(Hhat(:)))
        error('compute_main_nmse_2024_v2:InvalidEstimate', ...
            'Hhat must be finite (M_tau*G_nu*Nt)-by-Nr.');
    end
    required_channel = {'pathGains','pathDelays','pathDopplers','AoDs','AoAs'};
    for index = 1:numel(required_channel)
        if ~isfield(chanParams, required_channel{index})
            error('compute_main_nmse_2024_v2:MissingChannelParameter', ...
                'chanParams.%s is required.', required_channel{index});
        end
    end
    count = numel(chanParams.pathDelays);
    for index = 1:numel(required_channel)
        values = chanParams.(required_channel{index});
        if numel(values) ~= count || any(~isfinite(values(:)))
            error('compute_main_nmse_2024_v2:InvalidChannel', ...
                'Physical path arrays must be finite and equally sized.');
        end
    end
    if any(chanParams.pathDelays(:) < 0) || ...
            any(chanParams.pathDelays(:) ~= round(chanParams.pathDelays(:))) || ...
            max(chanParams.pathDelays) > params.cp_length_2024
        error('compute_main_nmse_2024_v2:InvalidDelay', ...
            'Physical delays must be integer taps covered by the CP.');
    end
    if ~isnumeric(data_locations) || numel(data_locations) ~= params.N || ...
            any(data_locations(:) < 1) || ...
            any(data_locations(:) > params.total_columns_2024) || ...
            any(diff(data_locations(:)) <= 0)
        error('compute_main_nmse_2024_v2:InvalidDataLocations', ...
            'data_locations must contain N increasing augmented-frame slots.');
    end
end
