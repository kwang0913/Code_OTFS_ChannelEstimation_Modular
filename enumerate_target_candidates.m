function [targets, rows] = enumerate_target_candidates(params)
%ENUMERATE_TARGET_CANDIDATES Feasible range-leg pairs on params.delay_resolution grid.
%   Deterministic: no RNG or file access. Inclusive [lo, hi] bounds:
%   params.delay_range (total delay, samples), params.aoa_range and params.aod_range (degrees).
%   ROWS is [DoD, DoA, Sign, totalDelay, x, y, AoD, AoA]. Angles are degrees.
%   Candidates lie between endpoint x coordinates (inclusive), on either
%   side of the Tx/Rx baseline. Degenerate collinear paths are excluded;
%   the direct LoS path is added separately by channel generation.
    check_range(params.delay_range,'delay_range',[0 Inf]);
    check_range(params.aoa_range,'aoa_range',[-90 90]);
    check_range(params.aod_range,'aod_range',[-90 90]);
    fsamp = params.M*params.df;
    tx = [params.Coord.Tx.x;params.Coord.Tx.y];
    rx = [params.Coord.Rx.x;params.Coord.Rx.y];
    x_min = min(tx(1),rx(1));
    x_max = max(tx(1),rx(1));
    R0 = norm(rx-tx);
    if R0 == 0
        error('enumerate_target_candidates:InvalidGeometry','Tx and Rx must be distinct.');
    end
    baseline_samples = R0*fsamp/params.c;
    geometry_tol = 16*eps(max(1,baseline_samples));
    baseline_axis = (rx-tx)/R0;
    basis = [baseline_axis,[-baseline_axis(2);baseline_axis(1)]];
    local.Tx = struct('x',0,'y',0);
    local.Rx = struct('x',R0,'y',0);
    ticks_per_sample = 1/params.delay_resolution;
    delay_ticks = round(ticks_per_sample*params.delay_range,12);
    min_total_tick = ceil(delay_ticks(1));
    max_total_tick = floor(delay_ticks(2));
    leg_ticks = max(0,max_total_tick-1);
    rows = zeros(10000,8); % Grow only for feasible candidates on the finer grid.
    count = 0;
    for dod_tick = 1:leg_ticks
        DoD = dod_tick/ticks_per_sample;
        for doa_tick = max(1,min_total_tick-dod_tick):(max_total_tick-dod_tick)
            DoA = doa_tick/ticks_per_sample;
            pathDelay = (dod_tick+doa_tick)/ticks_per_sample;
            % Strict triangle inequalities exclude direct/collinear paths.
            % Allow only roundoff-sized slack, not a physical separation gap.
            if pathDelay-baseline_samples <= geometry_tol || ...
                    baseline_samples-abs(DoD-DoA) <= geometry_tol
                continue;
            end
            [local_x,local_y] = compute_target_coord(DoD,DoA,local,fsamp,params.c);
            for Sign = [-1,1]
                point = tx + basis*[local_x;Sign*local_y];
                if point(1) < x_min || point(1) > x_max
                    continue;
                end
                aod = atan2d(point(2)-tx(2),abs(point(1)-tx(1)));
                aoa = atan2d(point(2)-rx(2),abs(point(1)-rx(1)));
                if aod < params.aod_range(1) || aod > params.aod_range(2) || ...
                        aoa < params.aoa_range(1) || aoa > params.aoa_range(2)
                    continue;
                end
                count = count+1;
                if count > size(rows,1)
                    rows(2*size(rows,1),8) = 0;
                end
                rows(count,:) = [DoD,DoA,Sign,pathDelay,point.',aod,aoa];
            end
        end
    end
    rows = sortrows(rows(1:count,:),[4,1,2,3]);
    targets = struct('DoD',num2cell(rows(:,1).'), ...
        'DoA',num2cell(rows(:,2).'),'Doppler',num2cell(zeros(1,count)), ...
        'Sign',num2cell(rows(:,3).'));
end

function check_range(value,name,limits)
    validateattributes(value,{'numeric'},{'real','finite','numel',2,'nondecreasing', ...
        '>=',limits(1),'<=',limits(2)},mfilename,name);
end
