function frame = generate_rg_frame_2022_v2(params)
%GENERATE_RG_FRAME_2022_V2 Build independent pilot and data RCP frames.
    required = {'M','N','Nt','ModOrder','N_pilots_2022', ...
        'data_power_2022','pilot_power_2022','total_columns_2022', ...
        'data_locations_2022','pilot_locations_2022','cp_length_2022'};
    for index = 1:numel(required)
        if ~isfield(params, required{index})
            error('generate_rg_frame_2022_v2:MissingParameter', ...
                'params.%s is required.', required{index});
        end
    end
    M = params.M; N = params.N; Nt = params.Nt;
    bits_per_symbol = log2(params.ModOrder);
    source_bits = randi([0 1], M*bits_per_symbol, N, Nt);
    symbol_indices = bit2int(reshape(source_bits, bits_per_symbol, []), ...
        bits_per_symbol);
    data_dd = sqrt(params.data_power_2022)*reshape( ...
        qammod(symbol_indices,params.ModOrder,UnitAveragePower=true), ...
        M,N,Nt);
    Fn = fft(eye(N))/sqrt(N);
    Fm = fft(eye(M))/sqrt(M);
    data_time = complex(zeros(M,N,Nt));
    pilot_tf = sqrt(params.pilot_power_2022)*qammod( ...
        randi([0 3],M,params.N_pilots_2022,Nt),4, ...
        UnitAveragePower=true);
    pilot_time = complex(zeros(size(pilot_tf)));
    for tx = 1:Nt
        data_time(:,:,tx) = data_dd(:,:,tx)*Fn';
        pilot_time(:,:,tx) = Fm'*pilot_tf(:,:,tx);
    end
    transmitted = complex(zeros(M,params.total_columns_2022,Nt));
    transmitted(:,params.pilot_locations_2022,:) = pilot_time;
    transmitted(:,params.data_locations_2022,:) = data_time;
    frame = struct('source_bits',source_bits,'data_dd',data_dd, ...
        'data_time',data_time,'pilot_tf',pilot_tf, ...
        'pilot_time',pilot_time,'pilots',pilot_time, ...
        'pilot_locations',params.pilot_locations_2022, ...
        'data_locations',params.data_locations_2022, ...
        'transmitted',transmitted, ...
        'pilot_rcp',local_add_rcp(pilot_time,params.cp_length_2022), ...
        'data_rcp',local_add_rcp(data_time,params.cp_length_2022));
end

function padded = local_add_rcp(samples,cp_length)
    [M,column_count,Nt] = size(samples);
    serial = reshape(samples,M*column_count,Nt);
    if cp_length == 0
        padded = serial;
    else
        padded = [serial(end-cp_length+1:end,:);serial];
    end
end
