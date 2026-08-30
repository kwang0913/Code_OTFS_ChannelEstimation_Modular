function padded = add_column_cp_2024_v2(payload, params)
%ADD_COLUMN_CP_2024_V2 Prepend an explicit cyclic prefix to every time column.
    if ~isnumeric(payload) || ndims(payload) > 3 || ...
            size(payload, 1) ~= params.M || ...
            size(payload, 2) ~= params.total_columns_2024 || ...
            size(payload, 3) ~= params.Nt || any(~isfinite(payload(:)))
        error('add_column_cp_2024_v2:InvalidPayload', ...
            'payload must be finite M-by-(N+Np)-by-Nt.');
    end
    cp = params.cp_length_2024;
    if cp > params.M
        error('add_column_cp_2024_v2:CPTooLong', ...
            'cp_length_2024 must not exceed M.');
    end
    if cp == 0
        padded = payload;
        return;
    end
    padded = cat(1, payload(end-cp+1:end, :, :), payload);
end
