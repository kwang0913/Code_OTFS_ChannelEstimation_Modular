function value = apply_ap_sip_dictionary_2023(descriptor, value, params, direction)
%APPLY_AP_SIP_DICTIONARY_2023 Apply a lazy AP-SIP sensing matrix or adjoint.
%   The descriptor represents the exact decoupled pilot/joint dictionary.
%   It retains only the affine waveform and decoder factors, never an
%   observation-by-channel-coefficient matrix.
    if ~isstruct(descriptor) || ~isfield(descriptor, 'kind') || ...
            ~strcmp(descriptor.kind, 'ap_sip_lazy')
        error('apply_ap_sip_dictionary_2023:InvalidDescriptor', ...
            'descriptor must be an ap_sip_lazy dictionary.');
    end
    if ~(ischar(direction) || (isstring(direction) && isscalar(direction)))
        error('apply_ap_sip_dictionary_2023:InvalidDirection', ...
            'direction must be forward or adjoint.');
    end
    direction = lower(char(direction));
    required = {'mode', 'XDD', 'row_count', 'column_count'};
    if ~all(isfield(descriptor, required)) || ...
            ~ismember(descriptor.mode, {'pilot', 'joint'})
        error('apply_ap_sip_dictionary_2023:InvalidDescriptor', ...
            'descriptor lacks the lazy sensing data.');
    end
    Q = params.M_tau * params.G_nu;
    if descriptor.column_count ~= Q * params.Nt || ...
            ~isequal(size(descriptor.XDD), [params.M, params.N, params.Nt])
        error('apply_ap_sip_dictionary_2023:InvalidDescriptor', ...
            'descriptor dimensions are incompatible with params.');
    end
    decoders = local_decoders(descriptor, params);
    if descriptor.row_count ~= params.M * sum(cellfun(@(x) size(x, 2), decoders))
        error('apply_ap_sip_dictionary_2023:InvalidDescriptor', ...
            'descriptor row count is incompatible with its decoders.');
    end
    switch direction
        case 'forward'
            if ~isnumeric(value) || ~isequal(size(value), [Q * params.Nt, 1]) || ...
                    any(~isfinite(value(:)))
                error('apply_ap_sip_dictionary_2023:InvalidInput', ...
                    'forward input must be a finite channel coefficient vector.');
            end
            H = complex(zeros(Q * params.Nt, params.Nr));
            H(:, 1) = value;
            [~, response] = apply_sparse_otfs_channel_2023(descriptor.XDD, ...
                struct('H', H), params, 0);
            value = local_decouple(response(:, :, 1), decoders);
        case 'adjoint'
            if ~isnumeric(value) || ~isequal(size(value), [descriptor.row_count, 1]) || ...
                    any(~isfinite(value(:)))
                error('apply_ap_sip_dictionary_2023:InvalidInput', ...
                    'adjoint input must be a finite decoupled observation vector.');
            end
            residual = local_recouple(value, decoders, params);
            value = local_channel_adjoint(descriptor.XDD, residual, params);
        otherwise
            error('apply_ap_sip_dictionary_2023:InvalidDirection', ...
                'direction must be forward or adjoint.');
    end
end

function decoders = local_decoders(descriptor, params)
    switch descriptor.mode
        case 'pilot'
            if ~isfield(descriptor, 'P') || ...
                    ~isequal(size(descriptor.P), [params.N, params.K2])
                error('apply_ap_sip_dictionary_2023:InvalidDescriptor', ...
                    'pilot descriptor requires its P decoder.');
            end
            decoders = {descriptor.P};
        case 'joint'
            if ~isfield(descriptor, 'D') || ~isfield(descriptor, 'P') || ...
                    ~isequal(size(descriptor.D), [params.N, params.N - params.K2]) || ...
                    ~isequal(size(descriptor.P), [params.N, params.K2])
                error('apply_ap_sip_dictionary_2023:InvalidDescriptor', ...
                    'joint descriptor requires D and P decoders.');
            end
            decoders = {descriptor.D, descriptor.P};
    end
end

function vector = local_decouple(response, decoders)
    vector = complex(zeros(size(response, 1) * sum(cellfun(@(x) size(x, 2), decoders)), 1));
    offset = 0;
    for index = 1:numel(decoders)
        block = response * decoders{index};
        count = numel(block);
        vector(offset + (1:count)) = block(:);
        offset = offset + count;
    end
end

function response = local_recouple(vector, decoders, params)
    response = complex(zeros(params.M, params.N));
    offset = 0;
    for index = 1:numel(decoders)
        width = size(decoders{index}, 2);
        count = params.M * width;
        block = reshape(vector(offset + (1:count)), params.M, width);
        response = response + block * decoders{index}';
        offset = offset + count;
    end
end

function coefficient = local_channel_adjoint(XDD, residual, params)
% Analytic adjoint of every delay/Doppler atom; no coefficient basis loop.
    M = params.M;
    N = params.N;
    Q = params.M_tau * params.G_nu;
    samples = (0:M-1).';
    coefficient = complex(zeros(Q * params.Nt, 1));
    for tx = 1:params.Nt
        for delay_index = 1:params.M_tau
            delay = params.delay_grid_2023(delay_index);
            shifted = circshift(XDD(:, :, tx), delay, 1);
            sample_inner = sum(conj(shifted) .* residual, 2);
            phase = exp(1j * 2 * pi * (samples - delay) * ...
                params.doppler_grid_2023(:).' / (M * N));
            rows = (tx - 1) * Q + (0:params.G_nu-1) * params.M_tau + delay_index;
            coefficient(rows) = phase' * sample_inner;
        end
    end
end
