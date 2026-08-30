function [Hhat,state] = estimate_rgomp_2022_v2(Y,Psi,params,n0)
%ESTIMATE_RGOMP_2022_V2 Algorithm 1 simultaneous row-group OMP.
    [Q,Nt,Nr,max_groups,tolerance,noise_floor_factor] = ...
        local_validate(Y,Psi,params,n0);
    Hhat = complex(zeros(Q*Nt,Nr));
    residual = Y;
    previous_energy = norm(residual,'fro')^2;
    noise_energy = numel(Y)*n0;
    noise_floor = noise_floor_factor*noise_energy;
    active_groups = zeros(0,1);
    active_columns = zeros(0,1);
    improvements = zeros(min(max_groups,Q),1);
    iteration = 0;
    stopping_reason = 'max_groups';
    if previous_energy <= noise_floor
        stopping_reason = 'noise_floor';
    else
        for iteration = 1:min(max_groups,Q)
            correlation = Psi'*residual;
            group_energy = zeros(Q,1);
            for q = 1:Q
                rows = q+(0:Nt-1)*Q;
                group_energy(q) = sum(abs(correlation(rows,:)).^2,'all');
            end
            group_energy(active_groups) = -Inf;
            [~,group] = max(group_energy);
            active_groups(end+1,1) = group; %#ok<AGROW>
            active_columns = reshape(active_groups+(0:Nt-1)*Q,[],1);
            coefficients = Psi(:,active_columns)\Y;
            residual = Y-Psi(:,active_columns)*coefficients;
            energy = norm(residual,'fro')^2;
            improvements(iteration) = max(previous_energy-energy,0);
            Hhat(:) = 0;
            Hhat(active_columns,:) = coefficients;
            if energy <= noise_floor
                stopping_reason = 'noise_floor';
                break;
            end
            if improvements(iteration) <= tolerance
                stopping_reason = 'improvement';
                break;
            end
            previous_energy = energy;
        end
    end
    state = struct('active_groups',active_groups,'active_columns',active_columns, ...
        'residual',residual,'residual_energy',norm(residual,'fro')^2, ...
        'noise_energy',noise_energy,'noise_floor',noise_floor, ...
        'energy_improvements',improvements(1:iteration), ...
        'iterations',iteration,'stopping_reason',stopping_reason, ...
        'converged',~strcmp(stopping_reason,'max_groups'));
end

function [Q,Nt,Nr,max_groups,tolerance,noise_floor_factor] = ...
        local_validate(Y,Psi,params,n0)
    if ~isnumeric(Y) || ~ismatrix(Y) || isempty(Y) || any(~isfinite(Y(:)))
        error('estimate_rgomp_2022_v2:InvalidObservation','Y must be finite and nonempty.');
    end
    if ~isnumeric(Psi) || size(Psi,1) ~= size(Y,1) || any(~isfinite(Psi(:)))
        error('estimate_rgomp_2022_v2:InvalidDictionary','Psi must match Y rows.');
    end
    required = {'Nt','Nr','rgomp_tolerance_2022','rgomp_max_groups_2022', ...
        'rgomp_noise_floor_factor_2022'};
    if ~isstruct(params) || ~all(isfield(params,required))
        error('estimate_rgomp_2022_v2:InvalidParams','Missing RG-OMP parameters.');
    end
    Nt = params.Nt; Nr = params.Nr;
    if size(Y,2) ~= Nr || mod(size(Psi,2),Nt) ~= 0
        error('estimate_rgomp_2022_v2:InvalidDimensions','Nt/Nr do not match Psi/Y.');
    end
    Q = size(Psi,2)/Nt;
    max_groups = params.rgomp_max_groups_2022;
    tolerance = params.rgomp_tolerance_2022;
    noise_floor_factor = params.rgomp_noise_floor_factor_2022;
    if max_groups < 1 || max_groups ~= round(max_groups) || ...
            tolerance < 0 || ~isfinite(tolerance) || ...
            ~isscalar(noise_floor_factor) || ~isreal(noise_floor_factor) || ...
            ~isfinite(noise_floor_factor) || noise_floor_factor <= 0 || ...
            ~isscalar(n0) || ~isreal(n0) || ~isfinite(n0) || n0 < 0
        error('estimate_rgomp_2022_v2:InvalidParams','Invalid stopping parameters.');
    end
end
