function test_2022_pipeline()
%TEST_2022_PIPELINE Exercise both estimators through main-physical detection.
    folder = fileparts(fileparts(mfilename('fullpath')));
    addpath(fileparts(folder),folder);
    target = struct('DoD',1,'DoA',1,'Doppler',0,'Sign',1);
    overrides = struct('M',8,'N',4,'Nt',2,'Nr',2,'ModOrder',4, ...
        'targets',target,'n_targets',[],'doppler_range',[], ...
        'skip_pilot_config',true,'N_pilots_2022',2,'M_tau',3,'G_nu',1, ...
        'doppler_oversampling_2022',1,'SNR_db_list',25,'trials',1, ...
        'em_max_iter_2022',1,'rgomp_max_groups_2022',3, ...
        'use_parallel',false,'verbose_output',true);
    params = config_params_2022_v2(overrides);
    oracle_overrides = overrides;
    oracle_overrides.Nt = 1;
    oracle_params = config_params_2022_v2(oracle_overrides);
    rng(31);
    [oracle_channel,~] = generate_channel_params(oracle_params);
    oracle_frame = generate_rg_frame_2022_v2(oracle_params);
    oracle_truth = generate_grid_truth_2022_v2(oracle_params,oracle_channel);
    [oracle_received,~] = apply_main_rg_channel_2022_v2( ...
        oracle_frame,oracle_channel,oracle_params,0);
    [oracle_ber,~,oracle_soft] = detect_lmmse_2022_v2( ...
        oracle_received,oracle_truth.H,oracle_frame,oracle_params,0);
    assert(oracle_ber == 0);
    assert(norm(oracle_soft(:)-oracle_frame.data_dd(:)) / ...
        norm(oracle_frame.data_dd(:)) < 1e-8);

    solver_params = config_params_2022_v2(struct( ...
        'M',512,'N',32,'Nt',2,'Nr',4,'ModOrder',4, ...
        'N_pilots_2022',8,'doppler_oversampling_2022',10, ...
        'SNR_db_list',15,'trials',1,'use_parallel',false, ...
        'verbose_output',false,'lmmse_max_iter_2022',100, ...
        'lmmse_tol_2022',1e-7));
    rng(9);
    [solver_channel,~] = generate_channel_params(solver_params);
    solver_frame = generate_rg_frame_2022_v2(solver_params);
    solver_truth = generate_grid_truth_2022_v2( ...
        solver_params,solver_channel);
    solver_n0 = solver_params.data_power_2022*10^(-15/10);
    [solver_received,~] = apply_main_rg_channel_2022_v2( ...
        solver_frame,solver_channel,solver_params,solver_n0);
    [~,~,~,solver_diagnostics] = detect_lmmse_2022_v2( ...
        solver_received,solver_truth.H,solver_frame,solver_params,solver_n0);
    assert(solver_diagnostics.converged);
    assert(solver_diagnostics.iterations <= solver_params.lmmse_max_iter_2022);
    assert(solver_diagnostics.relative_residual <= ...
        solver_params.lmmse_tol_2022);
    [output,results] = evalc('run_sim_2022_v2(params)');
    scalar_fields = {'NMSE_RGBL','NMSE_RGOMP','coefficient_NMSE_RGBL', ...
        'coefficient_NMSE_RGOMP','BER_RGBL','BER_RGOMP','runtime_seconds'};
    for index = 1:numel(scalar_fields)
        value = results.(scalar_fields{index});
        assert(isequal(size(value),[1 1]) && isfinite(value));
    end
    assert(results.BER_RGBL >= 0 && results.BER_RGBL <= 1);
    assert(results.BER_RGOMP >= 0 && results.BER_RGOMP <= 1);
    trial = results.trials{1};
    assert(trial.rgbl_state.iterations == 1);
    assert(trial.rgomp_state.iterations <= 3);
    assert(numel(trial.true_parameters.pathDelays) > 0);
    assert(contains(output,'Progress: 1/1'));
    assert(contains(output,'True physical paths:'));
    assert(~contains(output,'Target Coordinates'));
    progress_params = params;
    progress_params.trials = 2;
    progress_params.verbose_output = false;
    progress_output = evalc('run_sim_2022_v2(progress_params)');
    assert(contains(progress_output,'Progress: 2/2'));
    assert(~contains(progress_output,'True physical paths:'));
end
