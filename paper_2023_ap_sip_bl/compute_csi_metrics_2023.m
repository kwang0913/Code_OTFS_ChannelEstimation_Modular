function metrics = compute_csi_metrics_2023(Htrue, Hhat, supportTrue, lambda, threshold, params)
%COMPUTE_CSI_METRICS_2023 CSI/support scores with paper-facing operator NMSE.
    if nargin < 6
        params = [];
    end
    local_validate_inputs(Htrue, Hhat, supportTrue, lambda, threshold);

    coefficient_error_energy = norm(Htrue - Hhat, 'fro')^2;
    coefficient_reference_energy = norm(Htrue, 'fro')^2;
    coefficient_nmse = coefficient_error_energy / max(coefficient_reference_energy, eps);
    nmse = coefficient_nmse;
    if ~isempty(params)
        truth_operator = reconstruct_dd_channel_2023(Htrue, params);
        error_operator = reconstruct_dd_channel_2023(Htrue - Hhat, params);
        nmse = error_operator.frobenius_squared / max(truth_operator.frobenius_squared, eps);
    end

    supportTrue = logical(supportTrue(:));
    lambda = real(lambda(:));
    maximumLambda = max(lambda);
    if maximumLambda == 0
        estimatedSupport = false(size(lambda));
    else
        estimatedSupport = lambda >= maximumLambda * threshold;
    end
    truePositive = nnz(estimatedSupport & supportTrue);
    estimatedCount = nnz(estimatedSupport);
    trueCount = nnz(supportTrue);

    if estimatedCount == 0
        precision = 0;
    else
        precision = truePositive / estimatedCount;
    end
    if trueCount == 0
        recall = 0;
    else
        recall = truePositive / trueCount;
    end

    metrics = struct( ...
        'nmse', nmse, 'coefficient_nmse', coefficient_nmse, ...
        'estimated_support', logical(estimatedSupport), ...
        'support_precision', precision, ...
        'support_recall', recall);
end

function local_validate_inputs(Htrue, Hhat, supportTrue, lambda, threshold)
    if ~isnumeric(Htrue) || ~isnumeric(Hhat) || ~isequal(size(Htrue), size(Hhat)) || ...
            isempty(Htrue) || any(~isfinite(Htrue(:))) || any(~isfinite(Hhat(:)))
        error('compute_csi_metrics_2023:InvalidCSI', ...
            'Htrue and Hhat must be nonempty, finite numeric arrays of equal size.');
    end
    if ~(isnumeric(supportTrue) || islogical(supportTrue)) || ~isvector(supportTrue) || ...
            numel(supportTrue) ~= numel(lambda) || any(~isfinite(supportTrue(:)))
        error('compute_csi_metrics_2023:InvalidSupport', ...
            'supportTrue must be a finite vector with one entry per lambda value.');
    end
    if ~isnumeric(lambda) || ~isvector(lambda) || isempty(lambda) || ...
            ~isreal(lambda) || any(~isfinite(lambda(:))) || any(lambda(:) < 0)
        error('compute_csi_metrics_2023:InvalidLambda', ...
            'lambda must be a nonempty finite nonnegative real vector.');
    end
    if ~isscalar(threshold) || ~isreal(threshold) || ~isfinite(threshold) || ...
            threshold < 0 || threshold > 1
        error('compute_csi_metrics_2023:InvalidThreshold', ...
            'threshold must be a finite real scalar in [0, 1].');
    end
end
