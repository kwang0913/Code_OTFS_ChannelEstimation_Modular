function [Hhat, state] = estimate_pa_bl_mimo(Yp, Omega, Rw, params)
%ESTIMATE_PA_BL_MIMO Measurement-space row-group PA-BL.
%   Uses Woodbury's M-dimensional covariance and only materializes Sigma for
%   small regression systems.
    [Q,Nt,Nr,max_iterations,tolerance,floor_value] = local_validate_inputs(Yp,Omega,Rw,params);
    Rw=(Rw+Rw')/2; [~,noise_status]=chol(Rw,'lower');
    if noise_status~=0, error('estimate_pa_bl_mimo:NoiseCovarianceNotPD','Rw must be Hermitian positive definite.'); end
    lambda=ones(Q,1); converged=false; square_root_posterior_count=0;
    for iteration=1:max_iterations
        posterior=local_posterior(Omega,Yp,Rw,lambda,Nt,floor_value,false,false);
        square_root_posterior_count=square_root_posterior_count+ ...
            startsWith(posterior.mode,'square_root_');
        lambda_new=local_update_lambda(posterior.mean,posterior.diagonal,Q,Nt,Nr,floor_value);
        if any(~isfinite(lambda_new)), error('estimate_pa_bl_mimo:NonFiniteState','PA-BL variance update produced non-finite values.'); end
        squared_change=norm(lambda_new-lambda)^2; lambda=lambda_new;
        if squared_change<=tolerance, converged=true; break; end
    end
    materialize=Q*Nt<=1024;
    posterior=local_posterior(Omega,Yp,Rw,lambda,Nt,floor_value,materialize,true);
    square_root_posterior_count=square_root_posterior_count+ ...
        startsWith(posterior.mode,'square_root_');
    Hhat=posterior.mean;
    if any(~isfinite(Hhat(:)))||any(~isfinite(posterior.diagonal))||any(~isfinite(lambda)), error('estimate_pa_bl_mimo:NonFiniteState','PA-BL posterior contains non-finite values.'); end
    state=struct('lambda',lambda,'covariance_diagonal',posterior.diagonal,'prior_diagonal',posterior.prior_diagonal,'covariance_factor',posterior.factor,'factored_covariance',struct('diagonal',posterior.prior_diagonal,'factor',posterior.factor),'posterior_mode',posterior.mode,'square_root_posterior_count',square_root_posterior_count,'iterations',iteration,'converged',converged);
    if materialize, state.Sigma=posterior.full; end
end
function posterior=local_posterior(Omega,Y,Rw,lambda,Nt,floor_value,materialize,need_factor)
    d=repmat(max(lambda(:),floor_value),Nt,1);
    B=Omega.*d.';
    if numel(d)<=size(Omega,1)
        posterior=local_square_root_posterior( ...
            Omega,Y,Rw,d,B,materialize,need_factor);
        return;
    end
    S=B*Omega'+Rw;
    S=(S+S')/2;
    [L,status]=chol(S,'lower');
    if status~=0
        posterior=local_square_root_posterior( ...
            Omega,Y,Rw,d,B,materialize,need_factor);
        return;
    end
    solveY=L'\(L\Y);
    solveB=L'\(L\B);
    diagonal=max(real(d-sum(conj(B).*solveB,1).'),0);
    factor=[];
    if need_factor
        factor=B'/L';
    end
    posterior=struct('mean',B'*solveY,'diagonal',diagonal, ...
        'prior_diagonal',d,'factor',factor,'mode','measurement_cholesky');
    if materialize
        posterior.full=diag(d)-B'*solveB;
        posterior.full=(posterior.full+posterior.full')/2;
        posterior.full(1:numel(d)+1:end)=diagonal;
    end
end

function posterior=local_square_root_posterior( ...
        Omega,Y,Rw,d,B,materialize,need_factor)
    [noise_factor,noise_status]=chol(Rw,'lower');
    if noise_status~=0
        error('estimate_pa_bl_mimo:NoiseCovarianceNotPD', ...
            'Rw must be Hermitian positive definite.');
    end
    coefficient_count=numel(d);
    measurement_count=size(Omega,1);
    if coefficient_count<=measurement_count
        whitened_dictionary=noise_factor\Omega;
        whitened_observation=noise_factor\Y;
        augmented=[whitened_dictionary;diag(1./sqrt(d))];
        [~,R]=qr(augmented,0);
        mean_value=R\(R'\(whitened_dictionary'*whitened_observation));
        inverse_R=R\eye(coefficient_count);
        full_covariance=inverse_R*inverse_R';
        diagonal=max(real(sum(abs(inverse_R).^2,2)),0);
        factor=[];
        if need_factor
            reduction=diag(d)-full_covariance;
            reduction=(reduction+reduction')/2;
            [vectors,eigenvalues]=eig(reduction,'vector');
            factor=vectors.*sqrt(max(real(eigenvalues),0)).';
        end
        posterior=struct('mean',mean_value,'diagonal',diagonal, ...
            'prior_diagonal',d,'factor',factor,'mode','square_root_information');
        if materialize
            posterior.full=full_covariance;
            posterior.full(1:coefficient_count+1:end)=diagonal;
        end
        return;
    end

    square_root_dictionary=Omega.*sqrt(d).';
    covariance_root=[square_root_dictionary,noise_factor];
    [~,upper_factor]=qr(covariance_root',0);
    measurement_factor=upper_factor';
    solveY=measurement_factor'\(measurement_factor\Y);
    solveB=measurement_factor'\(measurement_factor\B);
    diagonal=max(real(d-sum(conj(B).*solveB,1).'),0);
    factor=[];
    if need_factor
        factor=B'/measurement_factor';
    end
    posterior=struct('mean',B'*solveY,'diagonal',diagonal, ...
        'prior_diagonal',d,'factor',factor,'mode','square_root_measurement');
    if materialize
        posterior.full=diag(d)-B'*solveB;
        posterior.full=(posterior.full+posterior.full')/2;
        posterior.full(1:coefficient_count+1:end)=diagonal;
    end
end
function lambda=local_update_lambda(mean_value,diagonal,Q,Nt,Nr,floor_value)
    lambda=zeros(Q,1);
    for q=1:Q
        rows=q+(0:Nt-1)*Q;
        lambda(q)=max(sum(abs(mean_value(rows,:)).^2,'all')/(Nt*Nr)+sum(diagonal(rows))/Nt,floor_value);
    end
end
function [Q,Nt,Nr,max_iterations,tolerance,floor_value]=local_validate_inputs(Yp,Omega,Rw,params)
    if ~isnumeric(Yp)||~ismatrix(Yp)||isempty(Yp)||any(~isfinite(Yp(:))), error('estimate_pa_bl_mimo:InvalidObservation','Yp must be a finite nonempty numeric matrix.'); end
    if ~isnumeric(Omega)||~ismatrix(Omega)||isempty(Omega)||any(~isfinite(Omega(:)))||size(Omega,1)~=size(Yp,1), error('estimate_pa_bl_mimo:InvalidDictionary','Omega must be finite and have one row per observation.'); end
    if ~isnumeric(Rw)||~ismatrix(Rw)||any(~isfinite(Rw(:)))||~isequal(size(Rw),[size(Yp,1),size(Yp,1)]), error('estimate_pa_bl_mimo:InvalidNoiseCovariance','Rw must be a finite square covariance matching Yp.'); end
    if norm(Rw-Rw','fro')>100*eps(max(1,norm(Rw,'fro'))), error('estimate_pa_bl_mimo:InvalidNoiseCovariance','Rw must be Hermitian.'); end
    if ~isstruct(params)||~isscalar(params)||~all(isfield(params,{'Nt','Nr','em_max_iter_2023','em_tol_2023','variance_floor_2023'})), error('estimate_pa_bl_mimo:InvalidParams','params is missing required PA-BL fields.'); end
    Nt=local_positive_integer(params.Nt,'Nt'); Nr=local_positive_integer(params.Nr,'Nr'); max_iterations=local_positive_integer(params.em_max_iter_2023,'em_max_iter_2023'); tolerance=local_nonnegative_scalar(params.em_tol_2023,'em_tol_2023'); floor_value=local_positive_scalar(params.variance_floor_2023,'variance_floor_2023');
    if size(Yp,2)~=Nr||mod(size(Omega,2),Nt)~=0, error('estimate_pa_bl_mimo:InvalidParams','Nr must match Yp columns and Omega columns must divide by Nt.'); end
    Q=size(Omega,2)/Nt;
end
function value=local_positive_integer(value,name)
    if ~isscalar(value)||~isnumeric(value)||~isreal(value)||~isfinite(value)||value<1||value~=round(value), error('estimate_pa_bl_mimo:InvalidParams','%s must be a positive integer.',name); end
end
function value=local_nonnegative_scalar(value,name)
    if ~isscalar(value)||~isnumeric(value)||~isreal(value)||~isfinite(value)||value<0, error('estimate_pa_bl_mimo:InvalidParams','%s must be a finite nonnegative scalar.',name); end
end
function value=local_positive_scalar(value,name)
    if ~isscalar(value)||~isnumeric(value)||~isreal(value)||~isfinite(value)||value<=0, error('estimate_pa_bl_mimo:InvalidParams','%s must be a finite positive scalar.',name); end
end
