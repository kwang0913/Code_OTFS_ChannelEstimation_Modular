function [aod, valid] = aod_from_geometry(aoa, delay_tap, M, df, c, Coord)
%AOD_FROM_GEOMETRY Derive AoD from AoA and total delay between endpoint x limits.
%   AoA/AoD are world-axis angles folded with abs(dx), matching generation.
%   Delay broadcasts with AoA. Coord contains the actual Tx/Rx endpoints.
    total = delay_tap/(M*df)*c;
    dx = Coord.Tx.x-Coord.Rx.x;
    dy = Coord.Tx.y-Coord.Rx.y;
    baseline = hypot(dx,dy);
    % Within the endpoint x interval the arrival ray points toward the other
    % endpoint in x; its y direction is determined by the signed AoA.
    ux = sign(dx)*cosd(aoa);
    uy = sind(aoa);
    projection = dx*ux+dy*uy;
    denominator = 2*(total-projection);
    valid = baseline>0 & total>baseline & abs(aoa)<=90 & denominator>0;
    denominator(~valid) = 1; % Invalid atoms are masked by the caller.
    arrival = (total.^2-baseline^2)./denominator;
    departure = total-arrival;
    x = Coord.Rx.x+arrival.*ux;
    y = Coord.Rx.y+arrival.*uy;
    valid = valid & arrival>0 & departure>0 & ...
        x>=min(Coord.Tx.x,Coord.Rx.x) & x<=max(Coord.Tx.x,Coord.Rx.x);
    aod = atan2d(y-Coord.Tx.y,abs(x-Coord.Tx.x));
    aod(~valid) = 0;
end
