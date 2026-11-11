function pb = default_pilot_blocks(params)
%DEFAULT_PILOT_BLOCKS Default structured pilot blocks (the value of pilot_blocks_def).
%   Tx1: delay block of length pilot_block_len ending at floor(M/2), at symbol floor(N/2).
%   Tx2: Doppler block of length pilot_block_len ending at the last Doppler bin, at the
%   subcarrier after the delay block. Tx3..Nt get n_pilotspertx random pilots each (config_pilots).
    blen = params.pilot_block_len;
    m_end = floor(params.M / 2);
    n_mid = floor(params.N / 2);
    if params.Nt < 2
        error('default_pilot_blocks:insufficientTransmitters', 'The default layout requires Nt >= 2.');
    end
    if blen > m_end || blen > params.N
        error('default_pilot_blocks:blockTooLong', 'pilot_block_len must not exceed floor(M/2) and N.');
    end
    pb = struct('dim', {'delay', 'doppler'}, ...
        'm', {m_end-(blen-1):m_end, m_end+1}, ...
        'n', {n_mid, params.N-(blen-1):params.N}, ...
        'tx', {1, 2}, ...
        'ref_m', {m_end-(blen-1), m_end+1}, ...
        'ref_n', {n_mid, params.N-(blen-1)}, ...
        'ref_tx', {1, 2});
end
