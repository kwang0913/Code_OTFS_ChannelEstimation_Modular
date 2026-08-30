function Psi = build_rg_dictionary_2022_v2(pilot_time,pilot_locations,params)
%BUILD_RG_DICTIONARY_2022_V2 Build the paper-style RCP pilot dictionary.
% Columns are delay-fast within Doppler and grouped by transmitter.
    Np = params.N_pilots_2022;
    if ~isnumeric(pilot_time) || size(pilot_time,1) ~= params.M || ...
            size(pilot_time,2) ~= Np || size(pilot_time,3) ~= params.Nt || ...
            any(~isfinite(pilot_time(:)))
        error('build_rg_dictionary_2022_v2:InvalidPilots', ...
            'pilot_time must be finite M-by-Np-by-Nt data.');
    end
    if ~isequal(reshape(pilot_locations,1,[]),params.pilot_locations_2022)
        error('build_rg_dictionary_2022_v2:InvalidLocations', ...
            'pilot_locations must match the independent training frame.');
    end
    sample_count = params.M*Np;
    Q = params.M_tau*params.G_nu;
    output_samples = params.cp_length_2022+(0:sample_count-1).';
    phase = exp(1j*2*pi*output_samples* ...
        params.doppler_grid_2022(:).' / ...
        params.doppler_reference_samples_2022);
    Psi = complex(zeros(sample_count,Q*params.Nt));
    for tx = 1:params.Nt
        pilot_serial = reshape(pilot_time(:,:,tx),sample_count,1);
        delayed = complex(zeros(sample_count,params.M_tau));
        for delay_index = 1:params.M_tau
            delayed(:,delay_index) = circshift(pilot_serial, ...
                params.delay_grid_2022(delay_index));
        end
        atoms = reshape(delayed,sample_count,params.M_tau,1).* ...
            reshape(phase,sample_count,1,params.G_nu);
        columns = (tx-1)*Q+(1:Q);
        Psi(:,columns) = reshape(atoms,sample_count,Q);
    end
end
