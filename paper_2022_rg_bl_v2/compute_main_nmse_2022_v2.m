function nmse = compute_main_nmse_2022_v2(chanParams,Hhat,params,data_locations)
%COMPUTE_MAIN_NMSE_2022_V2 Operator NMSE over original OTFS data slots.
    if ~isequal(reshape(data_locations,1,[]),params.data_locations_2022)
        error('compute_main_nmse_2022_v2:InvalidLocations', ...
            'data_locations must identify the independent data frame.');
    end
    denominator = params.doppler_reference_samples_2022;
    Q = params.M_tau*params.G_nu;
    physical_delays = chanParams.pathDelays(:);
    all_delays = unique([params.delay_grid_2022(:);physical_delays]);
    true_energy = 0; error_energy = 0;
    for data_index = 1:params.N
        samples = params.pilot_serial_length_2022+ ...
            params.cp_length_2022+(data_index-1)*params.M+ ...
            (0:params.M-1).';
        grid_phase = exp(1j*2*pi*samples* ...
            params.doppler_grid_2022(:).'/denominator);
        for tx = 1:params.Nt
            for rx = 1:params.Nr
                for delay = all_delays.'
                    true_gain = complex(zeros(params.M,1));
                    for path = find(physical_delays==delay).'
                        coefficient = chanParams.pathGains(path)* ...
                            exp(1j*2*pi*(rx-1)*params.ddlr*sind(chanParams.AoAs(path)))* ...
                            exp(-1j*2*pi*(tx-1)*params.ddlt*sind(chanParams.AoDs(path)));
                        true_gain = true_gain+coefficient*exp(1j*2*pi*samples* ...
                            chanParams.pathDopplers(path)/denominator);
                    end
                    estimated_gain = complex(zeros(params.M,1));
                    delay_index = find(params.delay_grid_2022==delay,1);
                    if ~isempty(delay_index)
                        rows = (tx-1)*Q+delay_index:params.M_tau:tx*Q;
                        estimated_gain = grid_phase*Hhat(rows,rx);
                    end
                    true_energy = true_energy+sum(abs(true_gain).^2);
                    error_energy = error_energy+sum(abs(estimated_gain-true_gain).^2);
                end
            end
        end
    end
    nmse = real(error_energy)/max(real(true_energy),eps);
    if ~isfinite(nmse) || nmse < 0
        error('compute_main_nmse_2022_v2:InvalidMetric','NMSE must be finite and nonnegative.');
    end
end
