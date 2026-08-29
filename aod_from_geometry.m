function [aod, valid] = aod_from_geometry(aoa, delay_tap, M, df, c, R0)
%AOD_FROM_GEOMETRY Derive AoD from (AoA, delay) via bistatic law of cosines.
%   aoa       - array of AoA candidates (degrees)
%   delay_tap - scalar or array of delay in taps (broadcast with aoa)
%   M, df, c, R0 - system parameters
%   Returns aod (degrees, same size as aoa) and valid (logical mask).

    range = delay_tap / (M * df) * c;       % total bistatic path length (m)

    % Physical validity: total path must exceed baseline
    valid = range > R0;

    % RoA = Range of Arrival (target→Rx), RoD = Range of Departure (Tx→target)
    % Derived from bistatic law of cosines: range^2 = R0^2 + 2*RoA*(range - R0*cos(aoa))
    RoA = (range.^2 - R0^2) ./ (2 * (range - R0 * cosd(abs(aoa))));
    RoD = range - RoA;

    valid = valid & (RoA > 0) & (RoD > 0);

    % Law of cosines on Tx–target–Rx triangle to get AoD magnitude
    cos_val = (R0^2 + RoD.^2 - RoA.^2) ./ (2 * R0 * RoD);
    cos_val = max(-1, min(1, cos_val));
    aod = acosd(cos_val) .* sign(aoa);     % sign(aoa): AoD mirrors AoA side of baseline

    aod(~valid) = 0;
end
