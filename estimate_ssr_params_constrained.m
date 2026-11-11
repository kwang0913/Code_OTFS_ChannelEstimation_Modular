function final_ests = estimate_ssr_params_constrained(YTF, XTF, coarse_ests, params, n0, chanParams, derivedParams)
%ESTIMATE_SSR_PARAMS_CONSTRAINED Law-of-cosines constrained SSR: AoD derived from (AoA,delay) geometry.
%   Two-stage search per target: (1) delay-Doppler grid, (2) AoA grid (AoD follows from geometry).
%   Binary refinement stops at the estimator minimum steps (ssr_min_steps,
%   defaults angle/delay/doppler_resolution); returned AoA/AoD, delays and
%   Dopplers lie on those grids, and gains are fitted with the returned values. If params.include_los,
%   the known LoS (last path, is_los) keeps the scenario-grid delay of
%   known_los_path so truth and receiver agree on it; otherwise all coarse paths move.
%   Gains use the pilot model y = c*Phi*gains + noise (pilot_noise_model).

    if params.include_los
        reflected = find(~coarse_ests.is_los);
        los = known_los_path(params);
    else
        reflected = 1:numel(coarse_ests.AoAs);
        los = struct('AoAs', [], 'AoDs', [], 'pathDelays', [], 'pathDopplers', []);
    end
    Kr = numel(reflected);  % movable paths
    K = Kr + params.include_los;
    Nr = params.Nr;
    M  = params.M;
    % Waveform timing and pilot layouts are centralized in the dictionary.
    df = params.df;
    pilot_blocks = params.pilot_blocks;

    n_meas = compute_n_meas(pilot_blocks, Nr);
    y = build_measurement_vector(YTF, XTF, pilot_blocks, Nr, n0);
    [c, s2] = pilot_noise_model(XTF, params, n0);  % y = c*Phi*gains + noise(s2)

    % --- Config ---
    if isfield(params, 'ssr_n_cand')
        n_cand = params.ssr_n_cand;
    else
        n_cand = 2;
    end
    min_step = ssr_min_steps(params);

    % Clamp the configured steps to the minimum steps BEFORE counting
    % intervals, so a minimum step above the configured step does not widen
    % the search range beyond ssr.range.
    cur_step.a = max(params.ssr.length.a, min_step.a);
    cur_step.r = max(params.ssr.length.r, min_step.r);
    cur_step.v = max(params.ssr.length.v, min_step.v);

    n_half.a = round(params.ssr.range.a / cur_step.a);
    n_half.r = round(params.ssr.range.r / cur_step.r);
    n_half.v = round(params.ssr.range.v / cur_step.v);

    % Reflections are movable; the LoS, if any, is fixed from endpoint geometry.
    centers = struct('AoAs', [coarse_ests.AoAs(reflected), los.AoAs], ...
                     'AoDs', [coarse_ests.AoDs(reflected), los.AoDs], ...
                     'pathDelays', [round(coarse_ests.pathDelays(reflected)/min_step.r)*min_step.r, los.pathDelays], ...
                     'pathDopplers', [snap(coarse_ests.pathDopplers(reflected), min_step.v), los.pathDopplers], ...
                     'is_los', [false(1, Kr), true(1, params.include_los)]);
    for i = 1:Kr
        [aod_i, valid_i] = aod_from_geometry(centers.AoAs(i), centers.pathDelays(i), M, df, params.c, params.Coord);
        if valid_i, centers.AoDs(i) = aod_i; end
    end
    final_ests = centers;

    % --- Binary refinement: halve step sizes each level until all below min_step ---
    level = 0;
    while true
        level = level + 1;
        cur_range.a = n_half.a * cur_step.a;
        cur_range.r = n_half.r * cur_step.r;
        cur_range.v = n_half.v * cur_step.v;

        if params.verbose_output
            fprintf('  BR level %d: step=[%.3f %.3f %.3f] range=[%.3f %.3f %.3f]\n', ...
                level, cur_step.a, cur_step.r, cur_step.v, cur_range.a, cur_range.r, cur_range.v);
        end

        % CoSaMP-style: compute residual per target, find best candidates
        Phi_current = c * build_response_matrix(centers, params, derivedParams);
        beta_current = (Phi_current'*Phi_current + s2*eye(K)) \ (Phi_current'*y);

        max_cands = Kr * n_cand + params.include_los;
        cand_columns = zeros(n_meas, max_cands);
        cand_params = repmat(struct('a',0,'d',0,'r',0,'v',0,'target',0), 1, max_cands);
        n_cand_total = 0;
        if params.include_los  % keep LoS in every joint fit without searching or pruning its geometry
            n_cand_total = 1;
            cand_columns(:, 1) = Phi_current(:, K);
            cand_params(1) = struct('a', los.AoAs, 'd', los.AoDs, ...
                'r', los.pathDelays, 'v', los.pathDopplers, 'target', K);
        end

        for i = 1:Kr
            r_i = y - Phi_current * beta_current + Phi_current(:,i) * beta_current(i);  % residual excluding target i

            % Stage 1: search delay-Doppler grid with AoA fixed
            [Phi_rv, grids_rv] = build_rv_dict_constrained( ...
                centers.AoAs(i), centers.pathDelays(i), centers.pathDopplers(i), ...
                cur_range.r, cur_step.r, cur_range.v, cur_step.v, ...
                params, derivedParams);
            cn = vecnorm(Phi_rv);
            valid_rv = cn > 0;
            if ~any(valid_rv), continue; end
            score_rv = -inf(size(cn));
            score_rv(valid_rv) = abs(Phi_rv(:, valid_rv)' * r_i).' ./ cn(valid_rv);
            [~, best_rv] = max(score_rv);
            [ir, iv] = ind2sub([grids_rv.n_r, grids_rv.n_v], best_rv);
            r_star = grids_rv.r(ir);
            v_star = grids_rv.v(iv);

            % Stage 2: search AoA grid with delay-Doppler fixed; AoD derived geometrically
            [Phi_ad, grids_ad] = build_ad_dict_constrained( ...
                centers.AoAs(i), r_star, v_star, ...
                cur_range.a, cur_step.a, ...
                params, derivedParams);
            cn = vecnorm(Phi_ad);
            valid_ad = grids_ad.valid_a & cn > 0;
            corr_ad = -inf(size(cn));
            corr_ad(valid_ad) = abs(Phi_ad(:, valid_ad)' * r_i).' ./ cn(valid_ad);
            [~, sorted_idx] = sort(corr_ad, 'descend');

            n_pick = min(n_cand, sum(valid_ad));  % keep valid candidates only
            for cc = 1:n_pick
                ci = sorted_idx(cc);
                n_cand_total = n_cand_total + 1;
                cand_columns(:, n_cand_total) = c * Phi_ad(:, ci);
                cand_params(n_cand_total) = struct( ...
                    'a', grids_ad.a(ci), 'd', grids_ad.d(ci), ...
                    'r', r_star, 'v', v_star, 'target', i);
            end
        end

        % Trim + joint LS + prune
        cand_columns = cand_columns(:, 1:n_cand_total);
        cand_params  = cand_params(1:n_cand_total);

        if n_cand_total > 0
            x_ls = (cand_columns'*cand_columns + s2*eye(n_cand_total)) \ (cand_columns'*y);
            for i = 1:Kr
                idx = find([cand_params.target] == i);
                if isempty(idx), continue; end
                [~, best] = max(abs(x_ls(idx)));
                bi = idx(best);
                final_ests.AoAs(i)         = snap(cand_params(bi).a, min_step.a);
                final_ests.AoDs(i)         = snap(cand_params(bi).d, min_step.a);
                final_ests.pathDelays(i)   = cand_params(bi).r;
                final_ests.pathDopplers(i) = snap(cand_params(bi).v, min_step.v);
            end
        end

        % Re-derive AoD, update centers
        for i = 1:Kr
            [aod_i, val_i] = aod_from_geometry(final_ests.AoAs(i), final_ests.pathDelays(i), M, df, params.c, params.Coord);
            if val_i, final_ests.AoDs(i) = snap(aod_i, min_step.a); end
        end
        centers = final_ests;

        if params.verbose_output
            fprintf('    -> AoA=[%s] delay=[%s] Doppler=[%s]\n', ...
                sprintf('%.1f ', final_ests.AoAs), ...
                sprintf('%.2f ', final_ests.pathDelays), ...
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
    Phi_est = c * build_response_matrix(final_ests, params, derivedParams);
    final_ests.pathGains = (Phi_est'*Phi_est + s2*eye(K)) \ (Phi_est'*y);
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

function x = snap(x, res)
% Round to the resolution grid (12-digit rounding removes multiplication noise).
    x = round(round(x / res) * res, 12);
end

function [Phi, grids] = build_rv_dict_constrained(center_a, center_r, center_v, ...
                                                   range_r, step_r, range_v, step_v, ...
                                                   params, derivedParams)
    lo = max(0, center_r - range_r);
    hi = min(params.M - eps(params.M), center_r + range_r);
    hi = min(hi, derivedParams.padLen);
    % Include grid points nearest integer symbol-overlap boundaries.
    resolution = ssr_min_steps(params).r;
    grid_r = unique(round([center_r-range_r:step_r:center_r+range_r, center_r, ceil(lo):floor(hi)]/resolution)*resolution);
    grid_r = grid_r(grid_r >= lo & grid_r <= hi);
    grid_v = center_v - range_v : step_v : center_v + range_v;
    n_r = numel(grid_r); n_v = numel(grid_v);
    [aod_arr, valid_r] = aod_from_geometry(center_a, grid_r, params.M, params.df, params.c, params.Coord);
    grids.r = grid_r; grids.v = grid_v; grids.n_r = n_r; grids.n_v = n_v;
    atoms = struct('AoAs', repmat(center_a, 1, n_r*n_v), ...
        'AoDs', repmat(aod_arr, 1, n_v), 'pathDelays', repmat(grid_r, 1, n_v), ...
        'pathDopplers', repelem(grid_v, n_r));
    Phi = build_response_matrix(atoms, params, derivedParams);
    Phi(:, ~repmat(valid_r, 1, n_v)) = 0;
end

function [Phi, grids] = build_ad_dict_constrained(center_a, fixed_r, fixed_v, ...
                                                    range_a, step_a, ...
                                                    params, derivedParams)
    grid_a = center_a - range_a : step_a : center_a + range_a;
    n_a = numel(grid_a);
    [grid_d, valid_a] = aod_from_geometry(grid_a, fixed_r, params.M, params.df, params.c, params.Coord);
    valid_a = valid_a & abs(grid_a) <= 90;
    grids.a = grid_a; grids.d = grid_d; grids.n_a = n_a; grids.valid_a = valid_a;
    atoms = struct('AoAs', grid_a, 'AoDs', grid_d, ...
        'pathDelays', repmat(fixed_r, 1, n_a), 'pathDopplers', repmat(fixed_v, 1, n_a));
    Phi = build_response_matrix(atoms, params, derivedParams);
    Phi(:, ~valid_a) = 0;
end
