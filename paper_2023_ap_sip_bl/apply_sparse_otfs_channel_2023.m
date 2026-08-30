function [Ydd, clean] = apply_sparse_otfs_channel_2023(XDD, truth, params, n0)
%APPLY_SPARSE_OTFS_CHANNEL_2023 Apply sparse MIMO delay-Doppler CSI.
%   Each atom is a circular delay shift followed by Doppler modulation.  The
%   dictionary builder invokes this operator, so both paths use this atom.

    local_validate_inputs(XDD, truth, params, n0);

    M = params.M;
    N = params.N;
    Nr = params.Nr;
    Q = params.M_tau * params.G_nu;
    clean = complex(zeros(M, N, Nr));

    samples = (0:M-1).';
    for tx = 1:params.Nt
        coefficients = reshape(truth.H((tx-1)*Q+(1:Q), :), Q, Nr);
        for delay_index = 1:params.M_tau
            delay = params.delay_grid_2023(delay_index);
            doppler_rows = delay_index:params.M_tau:Q;
            phase = exp(1j * 2 * pi * (samples-delay) * ...
                params.doppler_grid_2023(:).' / (M*N));
            received_coefficients = phase * coefficients(doppler_rows, :);
            shifted = circshift(XDD(:, :, tx), delay, 1);
            for rx = 1:Nr
                clean(:, :, rx) = clean(:, :, rx) + ...
                    received_coefficients(:, rx) .* shifted;
            end
        end
    end

    if n0 == 0
        Ydd = clean;
    else
        Ydd = clean + sqrt(n0 / 2) * (randn(M, N, Nr) + 1j * randn(M, N, Nr));
    end
end

function atom = local_dd_atom(symbols, delay, doppler)
%LOCAL_DD_ATOM Paper DD atom: delay-dependent Delta_l times row shift.
    shifted = circshift(symbols, [delay, 0]);
    source_samples = (0:size(symbols, 1)-1).' - delay;
    row_phase = exp(1j * 2 * pi * doppler * source_samples / ...
        (size(symbols, 1) * size(symbols, 2)));
    atom = row_phase .* shifted;
end

function local_validate_inputs(XDD, truth, params, n0)
    required = {'M', 'N', 'Nt', 'Nr', 'M_tau', 'G_nu', ...
        'delay_grid_2023', 'doppler_grid_2023'};
    for i = 1:numel(required)
        if ~isfield(params, required{i})
            error('apply_sparse_otfs_channel_2023:MissingParameter', ...
                'params.%s is required.', required{i});
        end
    end
    local_validate_dimensions(params);
    if ~isnumeric(XDD) || ~isequal(size(XDD), [params.M, params.N, params.Nt])
        error('apply_sparse_otfs_channel_2023:InvalidFrameSize', ...
            'XDD must be numeric M-by-N-by-Nt.');
    end
    if any(~isfinite(XDD(:)))
        error('apply_sparse_otfs_channel_2023:InvalidFrameValues', ...
            'XDD must contain only finite values.');
    end
    if ~isstruct(truth) || ~isfield(truth, 'H')
        error('apply_sparse_otfs_channel_2023:MissingChannel', ...
            'truth.H is required.');
    end
    Q = params.M_tau * params.G_nu;
    if ~isnumeric(truth.H) || ~isequal(size(truth.H), [Q * params.Nt, params.Nr])
        error('apply_sparse_otfs_channel_2023:InvalidChannelSize', ...
            'truth.H must be (M_tau*G_nu*Nt)-by-Nr.');
    end
    if any(~isfinite(truth.H(:)))
        error('apply_sparse_otfs_channel_2023:InvalidChannelValues', ...
            'truth.H must contain only finite values.');
    end
    local_validate_grids(params);
    if ~isscalar(n0) || ~isreal(n0) || ~isfinite(n0) || n0 < 0
        error('apply_sparse_otfs_channel_2023:InvalidNoiseVariance', ...
            'n0 must be a finite nonnegative real scalar.');
    end
end

function local_validate_dimensions(params)
    fields = {'M', 'N', 'Nt', 'Nr', 'M_tau', 'G_nu'};
    for i = 1:numel(fields)
        value = params.(fields{i});
        if ~isscalar(value) || ~isreal(value) || ~isfinite(value) || ...
                value < 1 || value ~= round(value)
            error('apply_sparse_otfs_channel_2023:InvalidDimension', ...
                'params.%s must be a finite positive integer.', fields{i});
        end
    end
end

function local_validate_grids(params)
    delay_grid = params.delay_grid_2023;
    if ~isnumeric(delay_grid) || ~isreal(delay_grid) || ...
            numel(delay_grid) ~= params.M_tau || any(~isfinite(delay_grid(:))) || ...
            any(delay_grid(:) ~= round(delay_grid(:)))
        error('apply_sparse_otfs_channel_2023:InvalidDelayGrid', ...
            'delay_grid_2023 must be finite, real, integer, and M_tau-long.');
    end
    doppler_grid = params.doppler_grid_2023;
    if ~isnumeric(doppler_grid) || ~isreal(doppler_grid) || ...
            numel(doppler_grid) ~= params.G_nu || any(~isfinite(doppler_grid(:)))
        error('apply_sparse_otfs_channel_2023:InvalidDopplerGrid', ...
            'doppler_grid_2023 must be finite, real, and G_nu-long.');
    end
end
