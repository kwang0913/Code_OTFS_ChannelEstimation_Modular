function final_ests = estimate_ssr_params_constrained(YTF, XTF, coarse_ests, params, n0, chanParams)
%ESTIMATE_SSR_PARAMS_CONSTRAINED Law-of-cosines constrained SSR: AoD derived from (AoA,delay) geometry.
%   Two-stage search per target: (1) delay-Doppler grid, (2) AoA grid (AoD follows from geometry).
%   Binary refinement halves step sizes until convergence, then final LMMSE for gains.

    K  = length(coarse_ests.AoAs);
    Nr = params.Nr;
    M  = params.M;
    N  = params.N;
    ddlr = params.ddlr;
    ddlt = params.ddlt;
    df = params.df;
    c  = params.c;
    R0 = params.init.R0;
    pilot_blocks = params.pilot_blocks;

    n_meas = compute_n_meas(pilot_blocks, Nr);
    y = build_measurement_vector(YTF, XTF, pilot_blocks, Nr);

    % --- Config ---
    if isfield(params, 'ssr_n_cand')
        n_cand = params.ssr_n_cand;
    else
        n_cand = 2;
    end
    min_step.a = 0.1;
    min_step.r = 0.2;
    min_step.v = 0.1;

    n_half.a = round(params.ssr.range.a / params.ssr.length.a);
    n_half.r = round(params.ssr.range.r / params.ssr.length.r);
    n_half.v = round(params.ssr.range.v / params.ssr.length.v);

    cur_step.a = params.ssr.length.a;
    cur_step.r = params.ssr.length.r;
    cur_step.v = params.ssr.length.v;

    % --- Initialize centers ---
    init_aod = zeros(1, K);
    for i = 1:K
        [aod_i, val_i] = aod_from_geometry(coarse_ests.AoAs(i), coarse_ests.pathDelays(i), M, df, c, R0);
        if val_i, init_aod(i) = aod_i;
        else,     init_aod(i) = coarse_ests.AoDs(i);
        end
    end

    centers = struct('AoAs', coarse_ests.AoAs, 'AoDs', init_aod, ...
                     'pathDelays', coarse_ests.pathDelays, ...
                     'pathDopplers', coarse_ests.pathDopplers);
    final_ests = struct('AoAs', zeros(1,K), 'AoDs', zeros(1,K), ...
                        'pathDelays', zeros(1,K), 'pathDopplers', zeros(1,K));

    % --- Binary refinement: halve step sizes each level until all below min_step ---
    level = 0;
    while true
        level = level + 1;
        cur_range.a = n_half.a * cur_step.a;
        cur_range.r = n_half.r * cur_step.r;
        cur_range.v = n_half.v * cur_step.v;

        if params.verbose_output
            fprintf('  BR level %d: step=[%.3f %.3f %.3f]\n', ...
                level, cur_step.a, cur_step.r, cur_step.v);
        end

        % CoSaMP-style: compute residual per target, find best candidates
        Phi_current = build_response_matrix(centers, params);
        beta_current = (Phi_current'*Phi_current + n0*eye(K)) \ (Phi_current'*y);

        max_cands = K * n_cand;
        cand_columns = zeros(n_meas, max_cands);
        cand_params = repmat(struct('a',0,'d',0,'r',0,'v',0,'target',0), 1, max_cands);
        n_cand_total = 0;

        for i = 1:K
            r_i = y - Phi_current * beta_current + Phi_current(:,i) * beta_current(i);  % residual excluding target i

            % Stage 1: search delay-Doppler grid with AoA fixed
            [Phi_rv, grids_rv] = build_rv_dict_constrained( ...
                centers.AoAs(i), centers.pathDelays(i), centers.pathDopplers(i), ...
                cur_range.r, cur_step.r, cur_range.v, cur_step.v, ...
                pilot_blocks, Nr, M, N, n_meas, ddlr, ddlt, df, c, R0);
            cn = vecnorm(Phi_rv); cn(cn==0) = Inf;
            [~, best_rv] = max(abs(Phi_rv' * r_i) ./ cn');
            [ir, iv] = ind2sub([grids_rv.n_r, grids_rv.n_v], best_rv);
            r_star = grids_rv.r(ir);
            v_star = grids_rv.v(iv);

            % Stage 2: search AoA grid with delay-Doppler fixed; AoD derived geometrically
            [Phi_ad, grids_ad] = build_ad_dict_constrained( ...
                centers.AoAs(i), r_star, v_star, ...
                cur_range.a, cur_step.a, ...
                pilot_blocks, Nr, M, N, n_meas, ddlr, ddlt, df, c, R0);
            cn = vecnorm(Phi_ad); cn(cn==0) = Inf;
            corr_ad = abs(Phi_ad' * r_i) ./ cn';
            [~, sorted_idx] = sort(corr_ad, 'descend');

            n_pick = min(n_cand, sum(grids_ad.valid_a));  % keep top-n_cand candidates per target
            for cc = 1:n_pick
                ci = sorted_idx(cc);
                n_cand_total = n_cand_total + 1;
                cand_columns(:, n_cand_total) = Phi_ad(:, ci);
                cand_params(n_cand_total) = struct( ...
                    'a', grids_ad.a(ci), 'd', grids_ad.d(ci), ...
                    'r', r_star, 'v', v_star, 'target', i);
            end
        end

        % Trim + joint LS + prune
        cand_columns = cand_columns(:, 1:n_cand_total);
        cand_params  = cand_params(1:n_cand_total);

        if n_cand_total > 0
            x_ls = (cand_columns'*cand_columns + n0*eye(n_cand_total)) \ (cand_columns'*y);
            for i = 1:K
                idx = find([cand_params.target] == i);
                if isempty(idx), continue; end
                [~, best] = max(abs(x_ls(idx)));
                bi = idx(best);
                final_ests.AoAs(i)         = round(cand_params(bi).a, 1);
                final_ests.AoDs(i)         = round(cand_params(bi).d, 1);
                final_ests.pathDelays(i)   = round(cand_params(bi).r);
                final_ests.pathDopplers(i) = cand_params(bi).v;
            end
        end

        % Re-derive AoD, update centers
        for i = 1:K
            [aod_i, val_i] = aod_from_geometry(final_ests.AoAs(i), final_ests.pathDelays(i), M, df, c, R0);
            if val_i, final_ests.AoDs(i) = round(aod_i, 1); end
        end
        centers = final_ests;

        if params.verbose_output
            fprintf('    -> AoA=[%s] delay=[%s] Doppler=[%s]\n', ...
                sprintf('%.1f ', final_ests.AoAs), ...
                sprintf('%.1f ', final_ests.pathDelays), ...
                sprintf('%.2f ', final_ests.pathDopplers));
        end

        % Termination
        if cur_step.a <= min_step.a && cur_step.r <= min_step.r && cur_step.v <= min_step.v
            break;
        end
        cur_step.a = max(cur_step.a / 2, min_step.a);
        cur_step.r = max(cur_step.r / 2, min_step.r);
        cur_step.v = max(cur_step.v / 2, min_step.v);
    end

    % --- Final LS ---
    for i = 1:K
        [aod_i, val_i] = aod_from_geometry(final_ests.AoAs(i), final_ests.pathDelays(i), M, df, c, R0);
        if val_i, final_ests.AoDs(i) = round(aod_i, 1); end
    end
    Phi_est = build_response_matrix(final_ests, params);
    final_ests.pathGains = (Phi_est'*Phi_est + n0*eye(K)) \ (Phi_est'*y);
    residual = y - Phi_est * final_ests.pathGains;
    final_ests.residual_nmse = real(residual' * residual) / real(y' * y);  % used by noise variance estimator as CSI error proxy

    if params.verbose_output
        fprintf('--- SSR Results ---\n');
        fprintf('AoA:     %s\n', sprintf('%5.2f ', final_ests.AoAs));
        fprintf('AoD:     %s\n', sprintf('%5.2f ', final_ests.AoDs));
        fprintf('Delay:   %s\n', sprintf('%5.2f ', final_ests.pathDelays));
        fprintf('Doppler: %s\n', sprintf('%5.2f ', final_ests.pathDopplers));
        fprintf('-------------------\n');
    end

    if params.show_figures
        gain_mag = abs(final_ests.pathGains(:)');
        labels = {'Delay [tap]', 'Doppler [tap]', 'AoA [°]', 'AoD [°]'};
        fields = {'pathDelays', 'pathDopplers', 'AoAs', 'AoDs'};
        true_fields = {'pathDelays', 'pathDopplers', 'AoAs', 'AoDs'};
        for fi = 1:4
            figure; hold on;
            stem(final_ests.(fields{fi}), gain_mag, 'filled', 'LineWidth', 1.5);
            for j = 1:params.targets_num
                xline(chanParams.(true_fields{fi})(j), 'r--', 'LineWidth', 2);
            end
            xlabel(labels{fi}, 'FontSize', 16); ylabel('Magnitude', 'FontSize', 16);
            set(gca, 'FontSize', 14, 'FontWeight', 'Bold'); grid on; hold off;
        end
    end
end


% =========================================================================
% Local functions
% =========================================================================

function [Phi, grids] = build_rv_dict_constrained(center_a, center_r, center_v, ...
                                                   range_r, step_r, range_v, step_v, ...
                                                   pilot_blocks, Nr, M, N, n_meas, ddlr, ddlt, df, c, R0)
    grid_r = center_r - range_r : step_r : center_r + range_r;
    grid_v = center_v - range_v : step_v : center_v + range_v;
    n_r = length(grid_r); n_v = length(grid_v);

    [aod_arr, valid_r] = aod_from_geometry(repmat(center_a, 1, n_r), grid_r, M, df, c, R0);
    grids.r = grid_r; grids.v = grid_v; grids.n_r = n_r; grids.n_v = n_v;

    a_fixed = exp(1j*2*pi*ddlr*(0:Nr-1).'*sind(center_a));
    k_arr = floor(grid_v); kappa_arr = grid_v - k_arr; nu_arr = grid_v;

    Phi = zeros(n_meas, n_r * n_v);
    idx = 0;
    for b = 1:length(pilot_blocks)
        blk = pilot_blocks(b);
        d_ref_arr   = exp(-1j*2*pi*ddlt*(blk.ref_tx-1)*sind(aod_arr));
        d_pilot_arr = exp(-1j*2*pi*ddlt*(blk.tx-1)*sind(aod_arr));

        n0_pos = blk.ref_n-1; m0_pos = blk.ref_m-1;
        xi_H = compute_xi_H(nu_arr, k_arr, kappa_arr, grid_r, n0_pos, m0_pos, N, M, 'grid');
        Phi(idx+(1:Nr), :) = a_fixed * reshape((xi_H .* d_ref_arr).', 1, []);
        idx = idx + Nr;

        for p = 1:blk.n_pilots
            if strcmp(blk.dim, 'random')
                n0_pos = blk.n(p)-1; m0_pos = blk.m(p)-1;
            elseif strcmp(blk.dim, 'delay')
                n0_pos = blk.n-1; m0_pos = blk.m(p)-1;
            else
                n0_pos = blk.n(p)-1; m0_pos = blk.m-1;
            end
            xi_H = compute_xi_H(nu_arr, k_arr, kappa_arr, grid_r, n0_pos, m0_pos, N, M, 'grid');
            Phi(idx+(1:Nr), :) = a_fixed * reshape((xi_H .* d_pilot_arr).', 1, []);
            idx = idx + Nr;
        end
    end

    for ir = find(~valid_r)
        Phi(:, ir + (0:n_v-1)*n_r) = 0;
    end
end


function [Phi, grids] = build_ad_dict_constrained(center_a, fixed_r, fixed_v, ...
                                                    range_a, step_a, ...
                                                    pilot_blocks, Nr, M, N, n_meas, ddlr, ddlt, df, c, R0)
    grid_a = center_a - range_a : step_a : center_a + range_a;
    n_a = length(grid_a);

    [grid_d, valid_a] = aod_from_geometry(grid_a, fixed_r, M, df, c, R0);
    grids.a = grid_a; grids.d = grid_d; grids.n_a = n_a; grids.valid_a = valid_a;

    component_a = exp(1j*2*pi*ddlr*(0:Nr-1).'*sind(grid_a));
    k_fixed = floor(fixed_v); kappa_fixed = fixed_v - k_fixed; nu_fixed = fixed_v;

    Phi = zeros(n_meas, n_a);
    idx = 0;
    for b = 1:length(pilot_blocks)
        blk = pilot_blocks(b);
        d_ref   = exp(-1j*2*pi*ddlt*(blk.ref_tx-1)*sind(grid_d));
        d_pilot = exp(-1j*2*pi*ddlt*(blk.tx-1)*sind(grid_d));

        n0_pos = blk.ref_n-1; m0_pos = blk.ref_m-1;
        xi_H = compute_xi_H(nu_fixed, k_fixed, kappa_fixed, fixed_r, n0_pos, m0_pos, N, M);
        Phi(idx+(1:Nr), :) = xi_H * (component_a .* d_ref);
        idx = idx + Nr;

        for p = 1:blk.n_pilots
            if strcmp(blk.dim, 'random')
                n0_pos = blk.n(p)-1; m0_pos = blk.m(p)-1;
            elseif strcmp(blk.dim, 'delay')
                n0_pos = blk.n-1; m0_pos = blk.m(p)-1;
            else
                n0_pos = blk.n(p)-1; m0_pos = blk.m-1;
            end
            xi_H = compute_xi_H(nu_fixed, k_fixed, kappa_fixed, fixed_r, n0_pos, m0_pos, N, M);
            Phi(idx+(1:Nr), :) = xi_H * (component_a .* d_pilot);
            idx = idx + Nr;
        end
    end
    Phi(:, ~valid_a) = 0;
end
