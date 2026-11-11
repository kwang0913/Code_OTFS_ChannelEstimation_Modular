function path = known_los_path(params)
%KNOWN_LOS_PATH Static direct-path geometry, without gain or phase information.
%   Angles use the root array convention: elevation is folded with abs(dx).
%   Delay uses the SCENARIO grid params.delay_resolution (never the estimator
%   step ssr.min_step.r): the LoS is known, not searched, so truth and
%   receiver must use the identical value or the NMSE would floor.
    dx = params.Coord.Rx.x - params.Coord.Tx.x;
    dy = params.Coord.Rx.y - params.Coord.Tx.y;
    delay = hypot(dx, dy) * params.M * params.df / params.c;
    path.pathDelays = round(delay/params.delay_resolution)*params.delay_resolution;
    path.pathDopplers = 0;
    path.AoDs = atan2d(dy, abs(dx));
    path.AoAs = atan2d(-dy, abs(dx));
    path.is_los = true;
end
