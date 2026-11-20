function trial = simulate_trial_2022_v2(params,SNR_db,trial_idx)
%SIMULATE_TRIAL_2022_V2 Run one main-physical RG-BL/RG-OMP trial.
%   Seeds rng(rng_seed_start + trial_idx - 1, 'threefry'). BER is NaN and detection
%   is skipped unless params.compute_ber (resolved by run_sim_2022_v2) is true.
    rng(params.rng_seed_start + trial_idx - 1, 'threefry');
    total_timer = tic;
    n0 = params.data_power_2022*10^(-SNR_db/10);
    chanParams = generate_channel_params(params);
    truth = generate_grid_truth_2022_v2(params,chanParams);
    frame = generate_rg_frame_2022_v2(params);
    [received,~] = apply_main_rg_channel_2022_v2(frame,chanParams,params,n0);
    Yp = reshape(received(:,frame.pilot_locations,:),[],params.Nr);
    Psi = build_rg_dictionary_2022_v2( ...
        frame.pilot_time,frame.pilot_locations,params);

    timer = tic;
    [Hrgbl,rgbl_state] = estimate_rgbl_2022_v2( ...
        Yp,Psi,max(n0,params.variance_floor_2022),params);
    rgbl_runtime_seconds = toc(timer);
    timer = tic;
    [Hrgomp,rgomp_state] = estimate_rgomp_2022_v2(Yp,Psi,params,n0);
    rgomp_runtime_seconds = toc(timer);

    NMSE_RGBL = compute_main_nmse_2022_v2( ...
        chanParams,Hrgbl,params,frame.data_locations);
    NMSE_RGOMP = compute_main_nmse_2022_v2( ...
        chanParams,Hrgomp,params,frame.data_locations);
    denominator = max(norm(truth.H,'fro')^2,eps);
    coefficient_NMSE_RGBL = norm(Hrgbl-truth.H,'fro')^2/denominator;
    coefficient_NMSE_RGOMP = norm(Hrgomp-truth.H,'fro')^2/denominator;
    support_RGBL = local_support_metrics(truth.support, ...
        rgbl_state.gamma >= params.support_threshold_2022*max(rgbl_state.gamma));
    estimated_rgomp = false(size(truth.support));
    estimated_rgomp(rgomp_state.active_groups) = true;
    support_RGOMP = local_support_metrics(truth.support,estimated_rgomp);

    if params.compute_ber
        timer = tic;
        [BER_RGBL,~,~,detection_rgbl] = detect_lmmse_2022_v2( ...
            received,Hrgbl,frame,params,n0);
        rgbl_detection_runtime = toc(timer);
        timer = tic;
        [BER_RGOMP,~,~,detection_rgomp] = detect_lmmse_2022_v2( ...
            received,Hrgomp,frame,params,n0);
        rgomp_detection_runtime = toc(timer);
    else
        BER_RGBL = NaN; BER_RGOMP = NaN;
        detection_rgbl = []; detection_rgomp = [];
        rgbl_detection_runtime = 0; rgomp_detection_runtime = 0;
    end
    [true_parameters,estimated_parameters] = local_parameter_summary( ...
        chanParams,truth,Hrgbl,Hrgomp,rgbl_state.gamma,params);

    trial = struct('NMSE_RGBL',NMSE_RGBL,'NMSE_RGOMP',NMSE_RGOMP, ...
        'coefficient_NMSE_RGBL',coefficient_NMSE_RGBL, ...
        'coefficient_NMSE_RGOMP',coefficient_NMSE_RGOMP, ...
        'BER_RGBL',BER_RGBL,'BER_RGOMP',BER_RGOMP, ...
        'support_RGBL',support_RGBL,'support_RGOMP',support_RGOMP, ...
        'H_RGBL',Hrgbl,'H_RGOMP',Hrgomp,'truth',truth, ...
        'rgbl_state',rgbl_state,'rgomp_state',rgomp_state, ...
        'detection_RGBL',detection_rgbl,'detection_RGOMP',detection_rgomp, ...
        'n0',n0,'rgbl_runtime_seconds',rgbl_runtime_seconds, ...
        'rgomp_runtime_seconds',rgomp_runtime_seconds, ...
        'detection_runtime_seconds',rgbl_detection_runtime+rgomp_detection_runtime, ...
        'runtime_seconds',toc(total_timer),'true_parameters',true_parameters, ...
        'estimated_parameters',estimated_parameters);
end

function metrics = local_support_metrics(truth,estimate)
    tp = sum(truth & estimate); fp = sum(~truth & estimate); fn = sum(truth & ~estimate);
    precision = tp/max(tp+fp,1); recall = tp/max(tp+fn,1);
    metrics = struct('precision',precision,'recall',recall, ...
        'f1',2*precision*recall/max(precision+recall,eps));
end

function [true_parameters,estimated] = local_parameter_summary( ...
        channel,truth,Hrgbl,Hrgomp,gamma,params)
    true_parameters = struct('pathDelays',channel.pathDelays(:).', ...
        'pathDopplers',channel.pathDopplers(:).','AoAs',channel.AoAs(:).', ...
        'AoDs',channel.AoDs(:).','pathGains',channel.pathGains(:).', ...
        'pathGridDelays',truth.pathGridDelays, ...
        'pathGridDopplers',truth.pathGridDopplers);
    estimated = struct('RGBL',local_rank(Hrgbl,gamma,params), ...
        'RGOMP',local_rank(Hrgomp,[],params));
end

function summary = local_rank(H,score,params)
    Q = params.M_tau*params.G_nu;
    norms = zeros(Q,1);
    for atom = 1:Q
        norms(atom) = norm(H(atom+(0:params.Nt-1)*Q,:),'fro');
    end
    if isempty(score), score = norms; end
    [~,ranked] = sort(score(:),'descend');
    count = min(params.targets_num,Q); atoms = ranked(1:count);
    delay_indices = mod(atoms-1,params.M_tau)+1;
    doppler_indices = floor((atoms-1)/params.M_tau)+1;
    summary = struct('atom_indices',atoms,'pathDelays', ...
        params.delay_grid_2022(delay_indices),'pathDopplers', ...
        params.doppler_grid_2022(doppler_indices),'score',score(atoms), ...
        'coefficient_norm',norms(atoms));
end
