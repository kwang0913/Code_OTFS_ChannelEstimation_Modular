% find_all_taps.m — Utility script: sweep (DoD,DoA,Sign) tuples to enumerate reachable targets.
% Maps each tuple to its resulting (AoD, AoA) via law-of-cosines geometry.
% Covers all integer AoD/AoA in [-60, 60].
close all;

params = config_params_std();
fsamp = params.M * params.df;
c = params.c;
Coord = params.Coord;
R0 = params.init.R0;  % Tx-Rx distance

max_tap = 200;

tap_map = containers.Map();

for DoD = 1:max_tap
    for DoA = 1:max_tap
        for Sign = [-1, 1]
            RoD = DoD * c / fsamp;
            RoA = DoA * c / fsamp;
            L = RoD + RoA;

            if (L < R0) || (abs(RoD - RoA) > R0)
                continue;
            end

            AoD_rad = acos((RoD^2 + R0^2 - RoA^2) / (2 * RoD * R0));
            tgt_x = Coord.Tx.x + RoD * cos(AoD_rad);
            tgt_y = (Coord.Tx.y + RoD * sin(AoD_rad)) * Sign;

            aod = round(atand(tgt_y / abs(tgt_x - Coord.Tx.x)),1);
            aoa = round(atand((tgt_y - Coord.Rx.y) / abs(tgt_x - Coord.Rx.x)),1);

            if aod >= -60 && aod <= 60 && aoa >= -60 && aoa <= 60
                key = sprintf('%d,%d', aod, aoa);  % key = "AoD,AoA" for grouping multiple (DoD,DoA) solutions
                entry = [DoD, DoA, Sign, DoD + DoA];
                if tap_map.isKey(key)
                    tap_map(key) = [tap_map(key); entry];
                else
                    tap_map(key) = entry;
                end
            end
        end
    end
end

% Load the packaged target pool independent of the caller's working directory.
project_root = fileparts(mfilename('fullpath'));
pool = jsondecode(fileread(fullfile(project_root, 'targets', 'pool_std.json')));
pool_set = containers.Map();
for pi = 1:numel(pool)
    pkey = sprintf('%d,%d,%d', pool(pi).DoD, pool(pi).DoA, pool(pi).Sign);
    pool_set(pkey) = true;
end

% Print summary (only entries with pathDelay < 12)
max_delay = 11;
keys = tap_map.keys();
fprintf('Reachable (AoD, AoA) pairs in [-60, 60] with pathDelay < %d:\n\n', max_delay);
for k = 1:length(keys)
    taps = tap_map(keys{k});
    taps = taps(taps(:,4) < max_delay, :);
    if isempty(taps)
        continue;
    end
    fprintf('AoD,AoA = (%s)  ->  %d solutions\n', keys{k}, size(taps, 1));
    for r = 1:size(taps, 1)
        DoD_r = taps(r,1); DoA_r = taps(r,2); Sign_r = taps(r,3);
        RoD_r = DoD_r * c / fsamp;
        RoA_r = DoA_r * c / fsamp;
        AoD_rad_r = acos((RoD_r^2 + R0^2 - RoA_r^2) / (2 * RoD_r * R0));
        x_r = Coord.Tx.x + RoD_r * cos(AoD_rad_r);
        y_r = (Coord.Tx.y + RoD_r * sin(AoD_rad_r)) * Sign_r;
        pkey = sprintf('%d,%d,%d', DoD_r, DoA_r, Sign_r);
        if pool_set.isKey(pkey)
            mark = ' [IN POOL]';
        else
            mark = '';
        end
        fprintf('    DoD=%3d, DoA=%3d, Sign=%+d, pathDelay=%d, x=%.2f, y=%.2f%s\n', ...
            DoD_r, DoA_r, Sign_r, taps(r,4), x_r, y_r, mark);
    end
end

%% Print copy-pasteable struct definition
fprintf('\n%% --- Copy-paste target definition ---\n');
fprintf('targets = struct( ...\n');
all_taps = [];
for k = 1:length(keys)
    taps = tap_map(keys{k});
    taps = taps(taps(:,4) < max_delay, :);
    if isempty(taps)
        continue;
    end
    [~, idx] = min(taps(:,4));
    all_taps = [all_taps; taps(idx,:)];
end
% Sort by pathDelay for readability
all_taps = sortrows(all_taps, 4);
n = size(all_taps, 1);
dod_str = strjoin(arrayfun(@(x) sprintf('%d', x), all_taps(:,1), 'UniformOutput', false), ', ');
doa_str = strjoin(arrayfun(@(x) sprintf('%d', x), all_taps(:,2), 'UniformOutput', false), ', ');
sgn_str = strjoin(arrayfun(@(x) sprintf('%+d', x), all_taps(:,3), 'UniformOutput', false), ', ');
fprintf('    ''DoD'',     {%s}, ...\n', dod_str);
fprintf('    ''DoA'',     {%s}, ...\n', doa_str);
fprintf('    ''Doppler'', {}, ...\n');
fprintf('    ''Sign'',    {%s}  ...\n', sgn_str);
fprintf(');\n');
fprintf('%% %d targets, pathDelay range: %d to %d\n\n', n, min(all_taps(:,4)), max(all_taps(:,4)));

%% Visualize scatterer positions on the x-y coordinate plane
scatter_x = [];  scatter_y = [];
delay_list = [];  aod_list = [];  aoa_list = [];
in_pool = [];

for k = 1:length(keys)
    taps = tap_map(keys{k});
    taps = taps(taps(:,4) < max_delay, :);
    if isempty(taps)
        continue;
    end
    [~, idx] = min(taps(:,4));
    DoD = taps(idx,1); DoA = taps(idx,2); Sign = taps(idx,3);
    RoD = DoD * c / fsamp;
    RoA = DoA * c / fsamp;
    AoD_rad = acos((RoD^2 + R0^2 - RoA^2) / (2 * RoD * R0));
    tgt_x = Coord.Tx.x + RoD * cos(AoD_rad);
    tgt_y = (Coord.Tx.y + RoD * sin(AoD_rad)) * Sign;
    scatter_x(end+1) = tgt_x;
    scatter_y(end+1) = tgt_y;
    delay_list(end+1) = taps(idx,4);
    angles = sscanf(keys{k}, '%f,%f');
    aod_list(end+1) = angles(1);
    aoa_list(end+1) = angles(2);
    pkey = sprintf('%d,%d,%d', DoD, DoA, Sign);
    in_pool(end+1) = pool_set.isKey(pkey);
end

figure; hold on;
% Not in pool: small gray circles
idx_out = ~in_pool;
scatter(scatter_x(idx_out), scatter_y(idx_out), 20, delay_list(idx_out), 'o');
% In pool: large filled circles with black edge
idx_in = logical(in_pool);
scatter(scatter_x(idx_in), scatter_y(idx_in), 80, delay_list(idx_in), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
for i = 1:length(scatter_x)
    label = sprintf(' (%g,%g)', aod_list(i), aoa_list(i));
    if in_pool(i)
        text(scatter_x(i), scatter_y(i), label, 'FontSize', 7, 'FontWeight', 'bold', 'Clipping', 'on');
    else
        text(scatter_x(i), scatter_y(i), label, 'FontSize', 6, 'Color', [.5 .5 .5], 'Clipping', 'on');
    end
end
plot(Coord.Tx.x, Coord.Tx.y, 'r^', 'MarkerSize', 12, 'LineWidth', 2);
plot(Coord.Rx.x, Coord.Rx.y, 'bs', 'MarkerSize', 12, 'LineWidth', 2);
colorbar;
colormap(jet);
xlabel('x (m)');
ylabel('y (m)');
title(sprintf('Scatterer positions (pathDelay < %d). Bold = in pool\\_std', max_delay));
legend('Not in pool', 'In pool', 'Tx', 'Rx', 'Location', 'best');
grid on;
axis equal;
