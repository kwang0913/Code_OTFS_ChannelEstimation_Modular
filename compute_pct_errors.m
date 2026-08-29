function [per_target_err, matched_ests] = compute_pct_errors(chanParams, ests)
%COMPUTE_PCT_ERRORS Per-target absolute errors of estimated vs true parameters.
%   [per_target_err, matched_ests] = compute_pct_errors(chanParams, ests)
%
%   per_target_err: Kx5 matrix, columns = [AoA, AoD, delay, doppler, gain].
%                   Rows aligned to true targets (matched by nearest AoA).
%                   Unmatched targets have NaN.
%   matched_ests:   Kx5 matrix of estimated values in the same target order.
%                   Unmatched targets have NaN.

    n_est  = length(ests.AoAs);
    n_true = length(chanParams.AoAs);
    n_match = min(n_est, n_true);

    per_target_err = NaN(n_true, 5);
    matched_ests   = NaN(n_true, 5);

    % Greedy match: assign each estimate to the closest unmatched true target by AoA
    available = 1:n_true;
    for ii = 1:n_match
        [~, best] = min(abs(chanParams.AoAs(available) - ests.AoAs(ii)));
        tidx = available(best);
        available(best) = [];
        per_target_err(tidx, :) = [ ...
            abs(chanParams.AoAs(tidx)        - ests.AoAs(ii)), ...
            abs(chanParams.AoDs(tidx)        - ests.AoDs(ii)), ...
            abs(chanParams.pathDelays(tidx)   - ests.pathDelays(ii)), ...
            abs(chanParams.pathDopplers(tidx) - ests.pathDopplers(ii)), ...
            abs(abs(chanParams.pathGains(tidx)) - abs(ests.pathGains(ii)))];
        matched_ests(tidx, :) = [ ...
            ests.AoAs(ii), ests.AoDs(ii), ...
            ests.pathDelays(ii), ests.pathDopplers(ii), ...
            abs(ests.pathGains(ii))];
    end
end
