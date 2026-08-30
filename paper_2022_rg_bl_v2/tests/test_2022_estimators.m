function test_2022_estimators()
%TEST_2022_ESTIMATORS Verify Algorithms 1 and 2 on controlled MMV systems.
    folder = fileparts(fileparts(mfilename('fullpath')));
    addpath(folder);
    rng(23);
    Q = 5; Nt = 2; Nr = 3; measurements = 16;
    [U,~] = qr(randn(measurements,Q*Nt)+1j*randn(measurements,Q*Nt),0);
    Psi = U;
    Htrue = complex(zeros(Q*Nt,Nr));
    active = [2 5];
    for q = active
        Htrue(q+(0:Nt-1)*Q,:) = randn(Nt,Nr)+1j*randn(Nt,Nr);
    end
    Y = Psi*Htrue;

    omp_params = struct('Nt',Nt,'Nr',Nr,'rgomp_tolerance_2022',1e-14, ...
        'rgomp_max_groups_2022',numel(active), ...
        'rgomp_noise_floor_factor_2022',1.1);
    [Homp,omp_state] = estimate_rgomp_2022_v2(Y,Psi,omp_params,0);
    assert(isequal(sort(omp_state.active_groups(:)),active(:)));
    assert(norm(Homp-Htrue,'fro')/norm(Htrue,'fro') < 1e-11);
    assert(omp_state.iterations == numel(active));

    noise_variance = 0.2;
    Rw = noise_variance*eye(measurements);
    Ynoisy = Y + sqrt(noise_variance/2)*(randn(size(Y))+1j*randn(size(Y)));
    noisy_omp_params = omp_params;
    noisy_omp_params.rgomp_max_groups_2022 = Q;
    noisy_omp_params.rgomp_noise_floor_factor_2022 = 1.25;
    [~,noisy_omp_state] = estimate_rgomp_2022_v2( ...
        Ynoisy,Psi,noisy_omp_params,noise_variance);
    assert(isequal(sort(noisy_omp_state.active_groups(:)),active(:)));
    assert(noisy_omp_state.iterations == numel(active));
    assert(strcmp(noisy_omp_state.stopping_reason,'noise_floor'));
    initial_covariance = inv(Psi'*(Rw\Psi)+eye(Q*Nt)); %#ok<MINV>
    initial_mean = initial_covariance*(Psi'*(Rw\Ynoisy));
    expected_gamma = zeros(Q,1);
    for q = 1:Q
        rows = q+(0:Nt-1)*Q;
        expected_gamma(q) = sum(abs(initial_mean(rows,:)).^2,'all')/(Nt*Nr) + ...
            sum(real(diag(initial_covariance(rows,rows))))/Nt;
    end
    one_step = struct('Nt',Nt,'Nr',Nr,'em_max_iter_2022',1, ...
        'em_tol_2022',0,'variance_floor_2022',1e-12);
    [~,one_step_state] = estimate_rgbl_2022_v2(Ynoisy,Psi,Rw,one_step);
    assert(norm(one_step_state.gamma-expected_gamma) < 1e-10);
    bl_params = struct('Nt',Nt,'Nr',Nr,'em_max_iter_2022',2, ...
        'em_tol_2022',0,'variance_floor_2022',1e-12);
    [Hbl,bl_state] = estimate_rgbl_2022_v2(Ynoisy,Psi,Rw,bl_params);
    assert(isequal(size(Hbl),[Q*Nt,Nr]));
    assert(isequal(size(bl_state.gamma),[Q,1]));
    assert(all(isfinite(Hbl(:))) && all(bl_state.gamma > 0));
    assert(isempty(bl_state.covariance_factor));
    d = repmat(bl_state.gamma,Nt,1);
    Sigma = inv(Psi'*(Rw\Psi)+diag(1./d)); %#ok<MINV>
    dense_mean = Sigma*(Psi'*(Rw\Ynoisy));
    assert(norm(Hbl-dense_mean,'fro')/max(norm(dense_mean,'fro'),eps) < 1e-10);
    assert(norm(diag(Sigma)-bl_state.covariance_diagonal) < 1e-10);
end
