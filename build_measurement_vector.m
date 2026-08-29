function y = build_measurement_vector(YTF, XTF, pilot_blocks, Nr)
%BUILD_MEASUREMENT_VECTOR Construct Y/X ratio vector from pilot blocks.
%   y = build_measurement_vector(YTF, XTF, pilot_blocks, Nr)
%
%   Layout: [ref(Nr) | pilot_1(Nr) | ... | pilot_P(Nr)] per block, concatenated.

    n_meas = compute_n_meas(pilot_blocks, Nr);
    y = zeros(n_meas, 1);
    idx = 0;
    for b = 1:length(pilot_blocks)
        blk = pilot_blocks(b);
        y(idx + (1:Nr)) = squeeze(YTF(blk.ref_m, blk.ref_n, 1:Nr)) ./ XTF(blk.ref_m, blk.ref_n, blk.ref_tx);
        idx = idx + Nr;
        for p = 1:blk.n_pilots
            if strcmp(blk.dim, 'random')
                y(idx + (1:Nr)) = squeeze(YTF(blk.m(p), blk.n(p), 1:Nr)) ./ XTF(blk.m(p), blk.n(p), blk.tx);
            elseif strcmp(blk.dim, 'delay')
                y(idx + (1:Nr)) = squeeze(YTF(blk.m(p), blk.n, 1:Nr)) ./ XTF(blk.m(p), blk.n, blk.tx);
            else
                y(idx + (1:Nr)) = squeeze(YTF(blk.m, blk.n(p), 1:Nr)) ./ XTF(blk.m, blk.n(p), blk.tx);
            end
            idx = idx + Nr;
        end
    end
end
