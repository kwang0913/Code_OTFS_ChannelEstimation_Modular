function test_2022_physical_model()
%TEST_2022_PHYSICAL_MODEL Verify independent paper-style RCP frames.
    folder = fileparts(fileparts(mfilename('fullpath')));
    project_root = fileparts(folder);
    addpath(project_root, folder);


    overrides = struct('M', 32, 'N', 4, 'Nt', 2, 'Nr', 2, ...
        'ModOrder', 4, 'N_pilots_2022', 2, 'M_tau', 12, 'G_nu', 25, ...
        'doppler_oversampling_2022', 2, 'SNR_db_list', 20, ...
        'trials', 1, 'use_parallel', false, 'verbose_output', false);
    params = config_params_2022_v2(overrides);
    assert(strcmp(params.padType, 'RCP'));
    assert(params.total_columns_2022 == 6);
    assert(isequal(params.pilot_locations_2022, 1:2));
    assert(isequal(params.data_locations_2022, 3:6));
    assert(isequal(params.delay_grid_2022(:), (0:11).'));
    assert(any(abs(params.doppler_grid_2022 - 0.5) < 10*eps));

    rng(17);
    frame = generate_rg_frame_2022_v2(params);
    assert(isequal(size(frame.pilot_tf), [32 2 2]));
    assert(isequal(size(frame.pilot_time), [32 2 2]));
    assert(isequal(size(frame.data_time), [32 4 2]));
    assert(isequal(size(frame.pilot_rcp), [64+params.cp_length_2022 2]));
    assert(isequal(size(frame.data_rcp), [128+params.cp_length_2022 2]));
    Fm = fft(eye(params.M))/sqrt(params.M);
    for tx = 1:params.Nt
        assert(norm(frame.pilot_time(:,:,tx)-Fm'*frame.pilot_tf(:,:,tx),'fro') ...
            < 1e-12);
        pilot_serial = reshape(frame.pilot_time(:,:,tx),[],1);
        data_serial = reshape(frame.data_time(:,:,tx),[],1);
        assert(norm(frame.pilot_rcp(:,tx)- ...
            [pilot_serial(end-params.cp_length_2022+1:end);pilot_serial]) < 1e-12);
        assert(norm(frame.data_rcp(:,tx)- ...
            [data_serial(end-params.cp_length_2022+1:end);data_serial]) < 1e-12);
    end

    channel = struct('pathGains', 0.7-0.2j, 'pathDelays', 1, ...
        'pathDopplers', 0.5, 'AoDs', 0, 'AoAs', 0);
    [received, clean, diagnostics] = ...
        apply_main_rg_channel_2022_v2(frame, channel, params, 0);
    assert(norm(received(:)-clean(:)) == 0);
    assert(diagnostics.pilot_serial_length == 64+params.cp_length_2022);
    assert(diagnostics.data_serial_length == 128+params.cp_length_2022);

    Yp = reshape(clean(:, frame.pilot_locations, :), [], params.Nr);
    Psi = build_rg_dictionary_2022_v2( ...
        frame.pilot_time, frame.pilot_locations, params);
    H = local_grid_coefficients(channel, params);
    relative_error = norm(Yp-Psi*H, 'fro') / max(norm(Yp, 'fro'), eps);
    assert(relative_error < 1e-11, ...
        'Physical pilot output and dictionary differ by %.3e.', relative_error);

    multipath = struct('pathGains', [0.7-0.2j, -0.3+0.1j], ...
        'pathDelays', [0, 2], 'pathDopplers', [0.5, -0.5], ...
        'AoDs', [0, 20], 'AoAs', [0, -15]);
    [multipath_received, multipath_clean] = ...
        apply_main_rg_channel_2022_v2(frame, multipath, params, 0);
    assert(isequal(size(multipath_received), ...
        [params.M, params.total_columns_2022, params.Nr]));
    multipath_Yp = reshape( ...
        multipath_clean(:,frame.pilot_locations,:),[],params.Nr);
    multipath_H = local_grid_coefficients(multipath,params);
    multipath_error = norm(multipath_Yp-Psi*multipath_H,'fro') / ...
        max(norm(multipath_Yp,'fro'),eps);
    assert(multipath_error < 1e-11, ...
        'Multipath physical output and dictionary differ by %.3e.', ...
        multipath_error);
end

function H = local_grid_coefficients(channel,params)
    Q = params.M_tau*params.G_nu;
    H = complex(zeros(Q*params.Nt,params.Nr));
    for path = 1:numel(channel.pathDelays)
        delay_index = find(params.delay_grid_2022==channel.pathDelays(path),1);
        doppler_index = find(abs(params.doppler_grid_2022- ...
            channel.pathDopplers(path)) < 10*eps,1);
        atom = delay_index+(doppler_index-1)*params.M_tau;
        ar = exp(1j*2*pi*(0:params.Nr-1).'*params.ddlr* ...
            sind(channel.AoAs(path)));
        at = exp(1j*2*pi*(0:params.Nt-1).'*params.ddlt* ...
            sind(channel.AoDs(path)));
        links = channel.pathGains(path)*ar*at';
        for tx = 1:params.Nt
            H((tx-1)*Q+atom,:) = H((tx-1)*Q+atom,:)+links(:,tx).';
        end
    end
end
