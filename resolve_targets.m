function params = resolve_targets(params, default_pool_file)
%RESOLVE_TARGETS Resolve the target list of a configuration (shared by std/HF).
%   params = RESOLVE_TARGETS(params, default_pool_file) fills target_pool,
%   target_pool_is_grid, selected_targets (when n_targets is given),
%   targets_num and the init.* fields from the override fields of PARAMS.
%
%   Priority: explicit targets > potential_targets > targets_file >
%   automatic grid (n_targets set) > packaged pool DEFAULT_POOL_FILE
%   (file name inside targets/). Selection by n_targets is deterministic
%   here (sizing only); the per-trial sample is drawn in generate_channel_params.

    root = fileparts(mfilename('fullpath'));
    params.target_pool_is_grid = false;
    if ~isempty(params.targets)
        target_list = params.targets;
    else
        if ~isempty(params.potential_targets)
            pool = params.potential_targets;
        elseif ~isempty(params.targets_file)
            targets_file = params.targets_file;
            if ~isfile(targets_file)
                targets_file = fullfile(root, targets_file);
            end
            if ~isfile(targets_file)
                error('resolve_targets:targetsFileNotFound', ...
                    'Target pool not found: %s', targets_file);
            end
            pool = jsondecode(fileread(targets_file));
        elseif ~isempty(params.n_targets)
            pool = enumerate_target_candidates(params);
            params.target_pool_is_grid = true;
        else
            pool = jsondecode(fileread(fullfile(root, 'targets', default_pool_file)));
        end
        params.target_pool = pool;
        if ~isempty(params.n_targets)
            validateattributes(params.n_targets,{'numeric'}, ...
                {'scalar','integer','nonnegative','<=',numel(pool)},mfilename,'n_targets');
            % Deterministic sizing only; actual selection is per-trial.
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
