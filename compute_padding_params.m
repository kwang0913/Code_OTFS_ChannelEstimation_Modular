function [Meff, numSamps, offset] = compute_padding_params(M, N, padLen, padType)
%COMPUTE_PADDING_PARAMS Compute Meff, numSamps, and offset from padding config.
%   [Meff, numSamps, offset] = compute_padding_params(M, N, padLen, padType)

    if strcmpi(padType, 'ZP') || strcmpi(padType, 'CP')
        Meff = M + padLen;  % ZP/CP extend each symbol by padLen samples
    else
        Meff = M;           % RCP/NONE: symbol length unchanged (prefix reuses tail)
    end

    if strcmpi(padType, 'RCP') || strcmpi(padType, 'CP')
        offset = padLen;    % receiver discards first padLen samples per symbol
    else
        offset = 0;
    end

    if strcmpi(padType, 'CP') || strcmpi(padType, 'ZP') || strcmpi(padType, 'NONE')
        numSamps = Meff * N;          % per-symbol padding: total = symbol_len * N
    else  % RCP, RZP
        numSamps = Meff * N + padLen;  % reduced CP: one shared prefix for entire frame
    end
end
