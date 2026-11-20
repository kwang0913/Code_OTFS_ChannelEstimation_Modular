% FIND_ALL_TAPS Enumerate independently sampled fractional range legs.
% DoD and DoA each use params.delay_resolution; no extra delay is added.
% Uses the same candidate enumeration and bounds as random trial sampling.
% LoS is reported separately, never exported as a reflection.
close all;
params = config_params_std();
Coord = params.Coord;
delay_range = params.delay_range;
[targets,rows] = enumerate_target_candidates(params);
count = size(rows,1);
% Only highlight a separately supplied pool, not the entire generated grid.
in_pool = false(count,1);
if isfield(params,'target_pool') && ~params.target_pool_is_grid
    pool = params.target_pool;
    pool_tuples = [[pool.DoD].', [pool.DoA].', [pool.Sign].'];
    in_pool = ismember(rows(:,1:3),pool_tuples,'rows');
end
all_taps = rows(:,1:4); % [DoD, DoA, Sign, total delay], including ALL solutions.
scatter_x = rows(:,5); scatter_y = rows(:,6);
aod_list = round(rows(:,7),1); aoa_list = round(rows(:,8),1);
delay_list = rows(:,4);
fprintf('All %d reflected tuples: AoD in [%.1f, %.1f] deg, AoA in [%.1f, %.1f] deg, pathDelay in [%.1f, %.1f] samples\n', ...
    count,params.aod_range,params.aoa_range,delay_range);
fprintf('Independent DoD/DoA grid step: %g sample\n\n',params.delay_resolution);
for k = 1:count
    mark = '';
    if in_pool(k), mark = ' [IN POOL]'; end
    fprintf('AoD=%5.1f, AoA=%5.1f: DoD=%.2f, DoA=%.2f, Sign=%+d, pathDelay=%.2f, x=%.2f, y=%.2f%s\n', ...
        aod_list(k),aoa_list(k),all_taps(k,1),all_taps(k,2),all_taps(k,3), ...
        delay_list(k),scatter_x(k),scatter_y(k),mark);
end
los = known_los_path(params);
fprintf('\nAutomatic LoS (not in targets): delay=%.2f, AoD=%.1f, AoA=%.1f, Doppler=0\n', ...
    los.pathDelays,los.AoDs,los.AoAs);

% Export all independently sampled range-leg pairs without extra offsets.
if count > 0
    dod_str = strjoin(compose('%.1f',all_taps(:,1)),', ');
    doa_str = strjoin(compose('%.1f',all_taps(:,2)),', ');
    sgn_str = strjoin(compose('%+d',all_taps(:,3)),', ');
    doppler_str = strjoin(compose('%d',zeros(count,1)),', ');
    fprintf('\n%% All candidates; select a subset before running a simulation.\n');
    fprintf('targets = struct( ...\n');
    fprintf('    ''DoD'',     {%s}, ...\n',dod_str);
    fprintf('    ''DoA'',     {%s}, ...\n',doa_str);
    fprintf('    ''Doppler'', {%s}, ...\n',doppler_str);
    fprintf('    ''Sign'',    {%s});\n',sgn_str);
    fprintf('%% Supply a subset as params.targets to preserve the enumerated range legs.\n');
else
    fprintf('\nNo reflected targets satisfy the configured bounds.\n');
end

if params.show_figures
    figure; hold on;
    h_out = scatter(scatter_x(~in_pool),scatter_y(~in_pool),20,delay_list(~in_pool),'o');
    h_in = scatter(scatter_x(in_pool),scatter_y(in_pool),80,delay_list(in_pool),'filled', ...
        'MarkerEdgeColor','k','LineWidth',1.5);
    % Data tips retain all tuple details without overlapping text labels.
    handles = [h_out,h_in]; groups = {~in_pool,in_pool};
    for group = 1:2
        mask = groups{group};
        handles(group).DataTipTemplate.DataTipRows = [ ...
            dataTipTextRow('DoD (samples)',all_taps(mask,1),'%.2f'), ...
            dataTipTextRow('DoA (samples)',all_taps(mask,2),'%.2f'), ...
            dataTipTextRow('Sign',all_taps(mask,3),'%+d'), ...
            dataTipTextRow('Delay (samples)',delay_list(mask),'%.2f'), ...
            dataTipTextRow('AoD (deg)',aod_list(mask),'%.1f'), ...
            dataTipTextRow('AoA (deg)',aoa_list(mask),'%.1f')];
    end
    h_tx = plot(Coord.Tx.x,Coord.Tx.y,'r^','MarkerSize',12,'LineWidth',2);
    h_rx = plot(Coord.Rx.x,Coord.Rx.y,'bs','MarkerSize',12,'LineWidth',2);
    colorbar; colormap(jet);
    xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('All reflected candidates (pathDelay in [%.1f, %.1f]). Filled = in pool',delay_range));
    legend([h_out,h_in,h_tx,h_rx],{'Other targets','In pool','Tx','Rx'},'Location','best');
    grid on; axis equal;
end
