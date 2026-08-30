function [Hhat,state] = estimate_rgbl_2022_v2(Y,Psi,Rw,params)
%ESTIMATE_RGBL_2022_V2 Algorithm 2 row-group Bayesian learning.
% A single gamma is shared by every transmitter and receiver at each atom.
    [Q,Nt,Nr,max_iterations,tolerance,floor_value] = ...
        local_validate(Y,Psi,Rw,params);
    [whitened_dictionary,whitened_observation] = local_whiten(Psi,Y,Rw);
    gram = whitened_dictionary'*whitened_dictionary;
    gram = (gram+gram')/2;
    rhs = whitened_dictionary'*whitened_observation;
    clear whitened_dictionary whitened_observation
    gamma = ones(Q,1);
    converged = false;
    squared_changes = zeros(max_iterations,1);
    for iteration = 1:max_iterations
        posterior = local_posterior(gram,rhs,gamma,Nt);
        gamma_new = local_gamma_update(posterior.mean,posterior.diagonal, ...
            Q,Nt,Nr,floor_value);
        squared_changes(iteration) = norm(gamma_new-gamma)^2;
        gamma = gamma_new;
        if squared_changes(iteration) <= tolerance
            converged = true;
            break;
        end
    end
    posterior = local_posterior(gram,rhs,gamma,Nt);
    Hhat = posterior.mean;
    if any(~isfinite(Hhat(:))) || any(~isfinite(gamma))
        error('estimate_rgbl_2022_v2:NonFiniteState', ...
            'RG-BL produced non-finite posterior state.');
    end
    state = struct('gamma',gamma, ...
        'covariance_diagonal',posterior.diagonal, ...
        'prior_diagonal',posterior.prior_diagonal, ...
        'covariance_factor',[], ...
        'iterations',iteration,'converged',converged, ...
        'squared_changes',squared_changes(1:iteration));
end

function posterior = local_posterior(gram,rhs,gamma,Nt)
    d = repmat(gamma(:),Nt,1);
    information = gram;
    diagonal_indices = 1:(size(information,1)+1):numel(information);
    information(diagonal_indices) = ...
        information(diagonal_indices)+reshape(1./d,1,[]);
    [R,status] = chol(information);
    if status ~= 0
        error('estimate_rgbl_2022_v2:PosteriorNotPD', ...
            'Posterior information matrix is not positive definite.');
    end
    mean_value = R\(R'\rhs);
    inverse_R = R\eye(numel(d));
    diagonal = max(real(sum(abs(inverse_R).^2,2)),0);
    posterior = struct('mean',mean_value,'diagonal',diagonal, ...
        'prior_diagonal',d);
end

function [WPsi,WY] = local_whiten(Psi,Y,Rw)
    if isscalar(Rw)
        scale = sqrt(Rw);
        WPsi = Psi/scale;
        WY = Y/scale;
    elseif isvector(Rw)
        scale = sqrt(Rw(:));
        WPsi = Psi./scale;
        WY = Y./scale;
    else
        covariance = (Rw+Rw')/2;
        [L,status] = chol(covariance,'lower');
        if status ~= 0
            error('estimate_rgbl_2022_v2:NoiseCovarianceNotPD', ...
                'Rw must be positive definite.');
        end
        WPsi = L\Psi;
        WY = L\Y;
    end
end

function gamma = local_gamma_update(mean_value,diagonal,Q,Nt,Nr,floor_value)
    gamma = zeros(Q,1);
    for q = 1:Q
        rows = q+(0:Nt-1)*Q;
        gamma(q) = max(sum(abs(mean_value(rows,:)).^2,'all')/(Nt*Nr) + ...
            sum(diagonal(rows))/Nt,floor_value);
    end
end

function [Q,Nt,Nr,max_iterations,tolerance,floor_value] = ...
        local_validate(Y,Psi,Rw,params)
    if ~isnumeric(Y) || ~ismatrix(Y) || isempty(Y) || any(~isfinite(Y(:)))
        error('estimate_rgbl_2022_v2:InvalidObservation','Y must be finite and nonempty.');
    end
    if ~isnumeric(Psi) || size(Psi,1) ~= size(Y,1) || any(~isfinite(Psi(:)))
        error('estimate_rgbl_2022_v2:InvalidDictionary','Psi must match Y rows.');
    end
    if ~isnumeric(Rw) || isempty(Rw) || any(~isfinite(Rw(:)))
        error('estimate_rgbl_2022_v2:InvalidNoise','Rw must be finite numeric data.');
    end
    if ~(isscalar(Rw) || (isvector(Rw) && numel(Rw)==size(Y,1)) || ...
            isequal(size(Rw),[size(Y,1),size(Y,1)]))
        error('estimate_rgbl_2022_v2:InvalidNoise','Rw dimensions do not match Y.');
    end
    if (isscalar(Rw) || isvector(Rw)) && ...
            (any(~isreal(Rw(:))) || any(Rw(:) <= 0))
        error('estimate_rgbl_2022_v2:InvalidNoise', ...
            'Scalar or diagonal noise variances must be positive real values.');
    end
    required = {'Nt','Nr','em_max_iter_2022','em_tol_2022','variance_floor_2022'};
    if ~isstruct(params) || ~all(isfield(params,required))
        error('estimate_rgbl_2022_v2:InvalidParams','Missing RG-BL parameters.');
    end
    Nt = params.Nt; Nr = params.Nr;
    if size(Y,2) ~= Nr || mod(size(Psi,2),Nt) ~= 0
        error('estimate_rgbl_2022_v2:InvalidDimensions','Nt/Nr do not match Psi/Y.');
    end
    Q = size(Psi,2)/Nt;
    max_iterations = params.em_max_iter_2022;
    tolerance = params.em_tol_2022;
    floor_value = params.variance_floor_2022;
    if max_iterations < 1 || max_iterations ~= round(max_iterations) || ...
            tolerance < 0 || floor_value <= 0 || ...
            any(~isfinite([max_iterations,tolerance,floor_value]))
        error('estimate_rgbl_2022_v2:InvalidParams','Invalid EM stopping parameters.');
    end
end
