function truth = generate_paper_channel_2023(params, chanParams)
%GENERATE_PAPER_CHANNEL_2023 Map geometric paths to sparse MIMO DD CSI.

    local_validate_params(params);
    path_delays = local_row_vector(chanParams, 'pathDelays');
    path_dopplers = local_row_vector(chanParams, 'pathDopplers');
    aoas = local_row_vector(chanParams, 'AoAs');
    aods = local_row_vector(chanParams, 'AoDs');
    path_gains = local_row_vector(chanParams, 'pathGains');
    K = numel(path_delays);
    if numel(path_dopplers) ~= K || numel(aoas) ~= K || ...
            numel(aods) ~= K || numel(path_gains) ~= K
        error('generate_paper_channel_2023:PathCountMismatch', ...
            'All physical path parameter arrays must have equal lengths.');
    end

    Q = params.M_tau * params.G_nu;
    truth.H = complex(zeros(Q * params.Nt, params.Nr));
    truth.support = false(Q, 1);
    truth.pathGains = path_gains;

    for k = 1:K
        q = local_atom_index(params, path_delays(k), path_dopplers(k));
        ar = exp(1j * 2 * pi * (0:params.Nr-1).' * params.ddlr * sind(aoas(k)));
        at = exp(1j * 2 * pi * (0:params.Nt-1).' * params.ddlt * sind(aods(k)));
        Hk = path_gains(k) * ar * at';

        for tx = 1:params.Nt
            rows = (tx - 1) * Q + q;
            truth.H(rows, :) = truth.H(rows, :) + Hk(:, tx).';
        end
        truth.support(q) = true;
    end
end

function local_validate_params(params)
    required = {'M_tau', 'G_nu', 'Nt', 'Nr', 'ddlt', 'ddlr', ...
        'delay_grid_2023', 'doppler_grid_2023'};
    for i = 1:numel(required)
        if ~isfield(params, required{i})
            error('generate_paper_channel_2023:MissingParameter', ...
                'params.%s is required.', required{i});
        end
    end
    if numel(params.delay_grid_2023) ~= params.M_tau || ...
            numel(params.doppler_grid_2023) ~= params.G_nu
        error('generate_paper_channel_2023:InvalidGrid', ...
            'Configured delay and Doppler grids must match M_tau and G_nu.');
    end
end

function values = local_row_vector(chanParams, field)
    if ~isfield(chanParams, field)
        error('generate_paper_channel_2023:MissingChannelParameter', ...
            'chanParams.%s is required.', field);
    end
    values = chanParams.(field)(:).';
    if any(~isfinite(values))
        error('generate_paper_channel_2023:InvalidChannelParameter', ...
            'chanParams.%s must be finite.', field);
    end
end


function q = local_atom_index(params, delay, doppler)
    if delay < params.delay_grid_2023(1) || delay > params.delay_grid_2023(end) || ...
            doppler < params.doppler_grid_2023(1) || doppler > params.doppler_grid_2023(end)
        error('generate_paper_channel_2023:OutOfGridPath', ...
            'Path delay and Doppler must lie within the configured dictionary grid.');
    end
    [~, delay_index] = min(abs(params.delay_grid_2023 - delay));
    [~, doppler_index] = min(abs(params.doppler_grid_2023 - doppler));
    q = (doppler_index - 1) * params.M_tau + delay_index;
end
