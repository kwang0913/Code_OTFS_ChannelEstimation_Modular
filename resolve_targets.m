function params = resolve_targets(params, default_pool_file)
%RESOLVE_TARGETS Resolve the configured target list and initial geometry.
%   PARAMS = RESOLVE_TARGETS(PARAMS, DEFAULT_POOL_FILE) applies the target
%   priority custom targets > random from pool > selected from pool and sets
%   target_pool (pool modes), selected_targets (n_targets mode), targets_num,
%   and init.{DoDs, DoAs, pathDelays, pathDopplers, Signs}.
%
%   The pool is potential_targets, else targets_file, else
%   targets/DEFAULT_POOL_FILE. A relative targets_file resolves against the
%   current folder first, then the package root.

    if ~isempty(params.targets)
        target_list = params.targets;
    else
        if ~isempty(params.potential_targets)
            pool = params.potential_targets;
        else
            package_root = fileparts(mfilename('fullpath'));
            if isempty(params.targets_file)
                targets_file = fullfile(package_root, 'targets', default_pool_file);
            else
                targets_file = params.targets_file;
                if ~isfile(targets_file)
                    targets_file = fullfile(package_root, targets_file);
                end
            end
            if ~isfile(targets_file)
                error('resolve_targets:targetsFileNotFound', ...
                    'Target pool not found: %s', targets_file);
            end
            pool = jsondecode(fileread(targets_file));
        end
        params.target_pool = pool;
        if ~isempty(params.n_targets)
            % Pool mode: use first n_targets as defaults; actual random
            % selection happens per-trial in generate_channel_params.
            params.selected_targets = 1:params.n_targets;
        end
        target_list = pool(params.selected_targets);
    end

    if ~isempty(params.target_dopplers)
        for i = 1:numel(target_list)
            target_list(i).Doppler = params.target_dopplers(i);
        end
    end

    params.targets_num = numel(target_list);
    params.init.DoDs         = [target_list.DoD];
    params.init.DoAs         = [target_list.DoA];
    params.init.pathDelays   = [target_list.DoD] + [target_list.DoA];  % bistatic delay = Tx-target + target-Rx path
    params.init.pathDopplers = [target_list.Doppler];
    params.init.Signs        = [target_list.Sign];
end
