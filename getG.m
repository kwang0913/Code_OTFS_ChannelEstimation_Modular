function G = getG(M, N, chanParams, padLen, padType)
%GETG Form time-domain channel matrix G from DD path parameters.

    [Meff, numSamps, ~] = compute_padding_params(M, N, padLen, padType);

    MN = Meff*N;
    P = length(chanParams.pathDelays);

    % Build per-delay-tap channel responses (complex sinusoids at Doppler freq)
    g = zeros(padLen+1, numSamps);
    for p = 1:P
        gp = chanParams.pathGains(p);
        lp = chanParams.pathDelays(p);
        vp = chanParams.pathDopplers(p);

        if strcmp(padType, 'CP') || strcmp(padType, 'RCP')
            % CP/RCP: phase origin offset by padLen to account for prefix samples
            g(lp+1,:) = g(lp+1,:) + gp*exp(1i*2*pi/MN * vp*((0:numSamps-1)-lp-padLen));
        else
            g(lp+1,:) = g(lp+1,:) + gp*exp(1i*2*pi/MN * vp*((0:numSamps-1)-lp));
        end
    end

    % Assemble G (sparse): collect (row, col, val) triplets per delay tap, build once
    unique_delays = unique(chanParams.pathDelays(:)).';
    all_rows = [];
    all_cols = [];
    all_vals = [];
    for l = unique_delays
        v = g(l+1, l+1:numSamps);
        rows = (l+1 : numSamps).';
        cols = (1 : numSamps-l).';
        all_rows = [all_rows; rows];
        all_cols = [all_cols; cols];
        all_vals = [all_vals; v(:)];
    end
    G = sparse(all_rows, all_cols, all_vals, numSamps, numSamps);
end
