function Omega = build_ap_sip_dictionary(symbols, precoder, params, mode)
%BUILD_AP_SIP_DICTIONARY Construct decoupled AP-SIP sensing matrices.
%   Pilot rows are vec(Y*P).  Joint rows stack vec(Y*D) then vec(Y*P).
%   Columns follow delay-fast, Doppler, transmitter truth.H ordering.

    mode = local_normalize_mode(mode);
    local_validate_params(params);
    switch mode
        case 'pilot'
            local_validate_numeric_contract(symbols, precoder, params, params.K2, ...
                'build_ap_sip_dictionary:InvalidPilotContract');
            XDD = local_affine_precoding(symbols, precoder);
            decoders = {precoder};
        case 'joint'
            [Xd, Xp, D, P] = local_validate_joint_contract(symbols, precoder, params);
            XDD = local_affine_precoding(Xd, D) + local_affine_precoding(Xp, P);
            decoders = {D, P};
    end

    Q = params.M_tau * params.G_nu;
    atom_count = Q * params.Nt;
    row_count = params.M * sum(cellfun(@(decoder) size(decoder, 2), decoders));
    if strcmp(mode, 'joint') && row_count * atom_count > 1e7
        % Only joint sensing can be lazy: PA-BL consumes a numeric pilot
        % matrix in its measurement-space posterior.
        Omega = struct('kind', 'ap_sip_lazy', 'row_count', row_count, ...
            'column_count', atom_count, 'mode', mode, 'XDD', XDD);
        if strcmp(mode, 'pilot')
            Omega.P = decoders{1};
        else
            Omega.D = decoders{1};
            Omega.P = decoders{2};
        end
        return;
    end
    Omega = local_analytic_dictionary(XDD, decoders, params, row_count, atom_count);
end

function Omega = local_analytic_dictionary(XDD, decoders, params, row_count, atom_count)
% Construct all atoms from shared delay shifts and Doppler phases.  This is
% algebraically identical to applying the channel to every coefficient basis
% vector, without invoking the full channel operator once per atom.
    M = params.M;
    Q = params.M_tau * params.G_nu;
    decoder_width = sum(cellfun(@(decoder) size(decoder, 2), decoders));
    samples = (0:M-1).';
    phase = exp(1j * 2 * pi * ...
        (reshape(samples, M, 1, 1, 1) - reshape(params.delay_grid_2023(:), 1, 1, params.M_tau, 1)) .* ...
        reshape(params.doppler_grid_2023(:), 1, 1, 1, params.G_nu) / (M * params.N));
    Omega = complex(zeros(row_count, atom_count));
    for tx = 1:params.Nt
        decoupled_shifts = complex(zeros(M, decoder_width, params.M_tau));
        for delay_index = 1:params.M_tau
            shifted = circshift(XDD(:, :, tx), params.delay_grid_2023(delay_index), 1);
            column_offset = 0;
            for decoder_index = 1:numel(decoders)
                width = size(decoders{decoder_index}, 2);
                columns = column_offset + (1:width);
                decoupled_shifts(:, columns, delay_index) = shifted * decoders{decoder_index};
                column_offset = column_offset + width;
            end
        end
        atoms = reshape(decoupled_shifts, M, decoder_width, params.M_tau, 1) .* phase;
        columns = (tx - 1) * Q + (1:Q);
        Omega(:, columns) = reshape(atoms, row_count, Q);
    end
end

function mode = local_normalize_mode(mode)
    if isstring(mode)
        mode = char(mode);
    end
    if ~ischar(mode) || ~ismember(lower(mode), {'pilot', 'joint'})
        error('build_ap_sip_dictionary:InvalidMode', ...
            'mode must be ''pilot'' or ''joint''.');
    end
    mode = lower(mode);
end

function [Xd, Xp, D, P] = local_validate_joint_contract(symbols, precoder, params)
    if ~isstruct(symbols) || ~isstruct(precoder) || ...
            ~isfield(symbols, 'Xd') || ~isfield(symbols, 'Xp') || ...
            ~isfield(precoder, 'D') || ~isfield(precoder, 'P')
        error('build_ap_sip_dictionary:InvalidJointContract', ...
            'joint mode requires symbols.Xd/Xp and precoder.D/P.');
    end
    Xd = symbols.Xd;
    Xp = symbols.Xp;
    D = precoder.D;
    P = precoder.P;
    local_validate_numeric_contract(Xd, D, params, params.N - params.K2, ...
        'build_ap_sip_dictionary:InvalidJointContract');
    local_validate_numeric_contract(Xp, P, params, params.K2, ...
        'build_ap_sip_dictionary:InvalidJointContract');
end

function local_validate_numeric_contract(symbols, precoder, params, expected_width, identifier)
    if ~isnumeric(symbols) || ndims(symbols) ~= 3 || ...
            size(symbols, 1) ~= params.M || size(symbols, 2) ~= expected_width || ...
            size(symbols, 3) ~= params.Nt || any(~isfinite(symbols(:)))
        error(identifier, 'symbols must be finite M-by-K-by-Nt data with the required K.');
    end
    if ~isnumeric(precoder) || ~ismatrix(precoder) || ...
            size(precoder, 1) ~= params.N || ...
            size(precoder, 2) ~= expected_width || any(~isfinite(precoder(:)))
        error(identifier, 'precoder must be finite N-by-K data with the required K.');
    end
end

function local_validate_params(params)
    required = {'M', 'N', 'Nt', 'Nr', 'M_tau', 'G_nu', 'K2', ...
        'delay_grid_2023', 'doppler_grid_2023'};
    for i = 1:numel(required)
        if ~isfield(params, required{i})
            error('build_ap_sip_dictionary:MissingParameter', ...
                'params.%s is required.', required{i});
        end
    end
    dimensions = {'M', 'N', 'Nt', 'Nr', 'M_tau', 'G_nu', 'K2'};
    for i = 1:numel(dimensions)
        value = params.(dimensions{i});
        if ~isscalar(value) || ~isreal(value) || ~isfinite(value) || ...
                value < 1 || value ~= round(value)
            error('build_ap_sip_dictionary:InvalidDimension', ...
                'params.%s must be a finite positive integer.', dimensions{i});
        end
    end
    if params.K2 >= params.N
        error('build_ap_sip_dictionary:InvalidDimension', ...
            'params.K2 must be smaller than params.N.');
    end
    delay_grid = params.delay_grid_2023;
    if ~isnumeric(delay_grid) || ~isreal(delay_grid) || ...
            numel(delay_grid) ~= params.M_tau || any(~isfinite(delay_grid(:))) || ...
            any(delay_grid(:) ~= round(delay_grid(:)))
        error('build_ap_sip_dictionary:InvalidDelayGrid', ...
            'delay_grid_2023 must be finite, real, integer, and M_tau-long.');
    end
    doppler_grid = params.doppler_grid_2023;
    if ~isnumeric(doppler_grid) || ~isreal(doppler_grid) || ...
            numel(doppler_grid) ~= params.G_nu || any(~isfinite(doppler_grid(:)))
        error('build_ap_sip_dictionary:InvalidDopplerGrid', ...
            'doppler_grid_2023 must be finite, real, and G_nu-long.');
    end
end

function XDD = local_affine_precoding(symbols, precoder)
    XDD = complex(zeros(size(symbols, 1), size(precoder, 1), size(symbols, 3)));
    for tx = 1:size(symbols, 3)
        XDD(:, :, tx) = symbols(:, :, tx) * precoder';
    end
end

function vector = local_decouple_response(response, decoders)
    row_count = size(response, 1) * sum(cellfun(@(decoder) size(decoder, 2), decoders));
    vector = complex(zeros(row_count, 1));
    offset = 0;
    for i = 1:numel(decoders)
        decoupled = response * decoders{i};
        count = numel(decoupled);
        vector(offset + 1:offset + count) = decoupled(:);
        offset = offset + count;
    end
end
