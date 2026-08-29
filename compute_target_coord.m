function [tgt_x, tgt_y] = compute_target_coord(DoD, DoA, Coord, fsamp, c)
%COMPUTE_TARGET_COORD Compute (x,y) target position from delay taps and geometry.

    RoD = DoD * c / fsamp;  % Range from Tx to target (meters)
    RoA = DoA * c / fsamp;  % Range from Rx to target (meters)
    L = RoD + RoA;

    Tx_x = Coord.Tx.x;  Tx_y = Coord.Tx.y;
    Rx_x = Coord.Rx.x;  Rx_y = Coord.Rx.y;
    R0 = sqrt((Rx_x - Tx_x)^2 + (Rx_y - Tx_y)^2);

    if (L < R0) || (abs(RoD - RoA) > R0)
        warning('Impossible geometry: circles do not intersect.');
        tgt_x = NaN;  tgt_y = NaN;
        return;
    end

    % Law of Cosines on triangle (Tx, Rx, Target)
    AoD_rad = acos((RoD^2 + R0^2 - RoA^2) / (2 * RoD * R0));

    tgt_x = Tx_x + RoD * cos(AoD_rad);
    tgt_y = Tx_y + RoD * sin(AoD_rad);
end
