function [received, clean, diagnostics] = ...
        apply_main_rg_channel_2022_v2(frame, chanParams, params, n0)
%APPLY_MAIN_RG_CHANNEL_2022_V2 Send pilot and data RCP frames continuously.
    required = {'pathGains','pathDelays','pathDopplers','AoDs','AoAs'};
    for index = 1:numel(required)
        if ~isfield(chanParams,required{index})
            error('apply_main_rg_channel_2022_v2:MissingChannelField', ...
                'chanParams.%s is required.',required{index});
        end
    end
    count = numel(chanParams.pathDelays);
    for index = 1:numel(required)
        values = chanParams.(required{index});
        if numel(values) ~= count || any(~isfinite(values(:)))
            error('apply_main_rg_channel_2022_v2:InvalidChannel', ...
                'Physical channel arrays must be finite and equally sized.');
        end
    end
    if any(chanParams.pathDelays(:) < 0) || ...
            any(chanParams.pathDelays(:) ~= round(chanParams.pathDelays(:))) || ...
            max(chanParams.pathDelays) > params.cp_length_2022
        error('apply_main_rg_channel_2022_v2:InvalidDelay', ...
            'Integer delays must be covered by cp_length_2022.');
    end
    if ~isscalar(n0) || ~isreal(n0) || ~isfinite(n0) || n0 < 0
        error('apply_main_rg_channel_2022_v2:InvalidNoise', ...
            'n0 must be finite and nonnegative.');
    end
    frame_required = {'pilot_rcp','data_rcp','pilot_locations','data_locations'};
    if ~isstruct(frame) || ~all(isfield(frame,frame_required))
        error('apply_main_rg_channel_2022_v2:InvalidFrame', ...
            'frame must contain independent pilot and data RCP blocks.');
    end
    pilot_serial_length = size(frame.pilot_rcp,1);
    data_serial_length = size(frame.data_rcp,1);
    if size(frame.pilot_rcp,2) ~= params.Nt || ...
            size(frame.data_rcp,2) ~= params.Nt
        error('apply_main_rg_channel_2022_v2:InvalidFrame', ...
            'RCP blocks must contain one column per transmitter.');
    end
    transmitted = [frame.pilot_rcp;frame.data_rcp];
    serial_length = pilot_serial_length+data_serial_length;
    fsamp = params.M*params.df;
    reference_duration = params.doppler_reference_samples_2022/fsamp;
    maximum_delay = max(chanParams.pathDelays);
    path_signals = complex(zeros(serial_length+maximum_delay,count));
    for path_index = 1:count
        path = struct('pathGains',chanParams.pathGains(path_index), ...
            'pathDelays',chanParams.pathDelays(path_index), ...
            'pathDopplerFreqs', ...
                chanParams.pathDopplers(path_index)/reference_duration);
        for tx = 1:params.Nt
            steering = exp(-1j*2*pi*(tx-1)*params.ddlt* ...
                sind(chanParams.AoDs(path_index)));
            atom = steering*dopplerChannel(transmitted(:,tx),fsamp,path);
            atom_length = numel(atom);
            path_signals(1:atom_length,path_index) = ...
                path_signals(1:atom_length,path_index)+atom;
        end
    end
    receive_steering = exp(1j*2*pi*(0:params.Nr-1).'*params.ddlr* ...
        sind(chanParams.AoAs(:).'));
    serial_full = (receive_steering*path_signals.').';
    serial_clean = serial_full(1:serial_length,:);
    cp = params.cp_length_2022;
    pilot_padded = serial_clean(1:pilot_serial_length,:);
    data_padded = serial_clean(pilot_serial_length+(1:data_serial_length),:);
    pilot_clean = reshape(pilot_padded(cp+1:end,:), ...
        params.M,params.N_pilots_2022,params.Nr);
    data_clean = reshape(data_padded(cp+1:end,:), ...
        params.M,params.N,params.Nr);
    clean = complex(zeros(params.M,params.total_columns_2022,params.Nr));
    clean(:,frame.pilot_locations,:) = pilot_clean;
    clean(:,frame.data_locations,:) = data_clean;
    received = clean+sqrt(n0/2)*(randn(size(clean))+1j*randn(size(clean)));
    diagnostics = struct('serial_clean_full',serial_full, ...
        'serial_clean',serial_clean,'pilot_padded_clean',pilot_padded, ...
        'data_padded_clean',data_padded,'fsamp',fsamp, ...
        'reference_duration',reference_duration, ...
        'pilot_serial_length',pilot_serial_length, ...
        'data_serial_length',data_serial_length);
end
