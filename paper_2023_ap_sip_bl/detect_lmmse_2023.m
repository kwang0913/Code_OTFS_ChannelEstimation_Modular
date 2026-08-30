function [XdHat,bitsHat] = detect_lmmse_2023(Ydd,Hhat,frame,params,n0,posterior_covariance)
%DETECT_LMMSE_2023 Decoupled AP-SIP per-column uncertainty-aware LMMSE.
% No frame-wide data action is formed: every data column solves one M*Nt
% channel block, whose factorization is reused across all N-K2 columns.
    if nargin<6, posterior_covariance=[]; end
    [M,N,Nt,Nr,data_columns,mod_order,data_power]=local_validate_inputs(Ydd,Hhat,frame,params,n0,posterior_covariance);
    decoupled=complex(zeros(M,data_columns,Nr));
    for rx=1:Nr, decoupled(:,:,rx)=Ydd(:,:,rx)*frame.D; end
    Hblock=local_channel_block(Hhat,params);
    xi=local_xi_block(posterior_covariance,params);
    regularization=max(n0/data_power,eps(max(1,norm(Hblock,'fro')^2)));
    output_covariance=(xi+xi')/2+regularization*eye(M);
    [output_factor,output_status]=chol(output_covariance,'lower');
    if output_status~=0, error('detect_lmmse_2023:OutputCovarianceNotPD','The uncertainty-aware output covariance must be positive definite.'); end
    normal=eye(M*Nt);
    for receiver=1:Nr
        rows=(receiver-1)*M+(1:M);
        receiver_block=Hblock(rows,:);
        normal=normal+receiver_block'*(output_factor'\(output_factor\receiver_block));
    end
    normal=(normal+normal')/2;
    [normal_factor,normal_status]=chol(normal,'lower');
    if normal_status~=0, error('detect_lmmse_2023:NormalMatrixNotPD','The uncertainty-aware data normal matrix must be positive definite.'); end
    XdHat=complex(zeros(M,data_columns,Nt));
    for column=1:data_columns
        rhs=complex(zeros(M*Nt,1));
        for receiver=1:Nr
            rows=(receiver-1)*M+(1:M);
            y=decoupled(:,column,receiver);
            rhs=rhs+Hblock(rows,:)'*(output_factor'\(output_factor\y));
        end
        x=normal_factor'\(normal_factor\rhs);
        XdHat(:,column,:)=reshape(x,M,1,Nt);
    end
    hard=qamdemod(XdHat(:)/sqrt(data_power),mod_order,UnitAveragePower=true);
    bitsHat=reshape(int2bit(hard,log2(mod_order)),size(frame.bits));
    if any(~isfinite(XdHat(:)))||any(~ismember(bitsHat(:),[0;1])), error('detect_lmmse_2023:NonFiniteOutput','LMMSE detection produced an invalid output.'); end
end
function Hblock=local_channel_block(H,params)
% Assemble the exact per-column channel block from cyclic delay diagonals.
    M=params.M; Nt=params.Nt; Nr=params.Nr; Q=params.M_tau*params.G_nu;
    rows_total=M*Nr; samples=(0:M-1).';
    Hblock=complex(zeros(rows_total,M*Nt));
    for tx=1:Nt
        column_offset=(tx-1)*M;
        for delay_index=1:params.M_tau
            delay=params.delay_grid_2023(delay_index);
            source=mod(samples-delay,M)+1;
            atom_rows=delay_index:params.M_tau:Q;
            phase=exp(1j*2*pi*(samples-delay)*params.doppler_grid_2023(:).'/(M*params.N));
            for rx=1:Nr
                values=phase*H((tx-1)*Q+atom_rows,rx);
                rows=(rx-1)*M+(1:M).';
                indices=rows+(column_offset+source-1)*rows_total;
                Hblock(indices)=Hblock(indices)+values;
            end
        end
    end
end
function xi=local_xi_block(covariance,params)
    QNt=params.M_tau*params.G_nu*params.Nt;
    if isempty(covariance), xi=complex(zeros(params.M)); return; end
    if isstruct(covariance)
        if ~isfield(covariance,'diagonal')||~isfield(covariance,'factor')
            error('detect_lmmse_2023:InvalidPosteriorCovariance','Factored covariance needs diagonal and factor.');
        end
        diagonal=covariance.diagonal;
        U=covariance.factor;
        if ~isequal(size(diagonal),[QNt,1])||size(U,1)~=QNt|| ...
                any(~isfinite(diagonal))||any(~isfinite(U(:)))||any(diagonal<0)
            error('detect_lmmse_2023:InvalidPosteriorCovariance','Invalid factored covariance.');
        end
        covariance=diag(diagonal)-U*U';
    elseif isscalar(covariance)
        covariance=covariance*eye(QNt);
    end
    if ~isnumeric(covariance)||~isequal(size(covariance),[QNt,QNt])||any(~isfinite(covariance(:)))
        error('detect_lmmse_2023:InvalidPosteriorCovariance','posterior covariance has incompatible dimensions.');
    end
    xi=local_output_covariance_xi((covariance+covariance')/2,params);
end
function xi=local_output_covariance_xi(covariance,params)
% Eq. (74): E{DeltaH*DeltaH^H} for one receive antenna.
    M=params.M; Nt=params.Nt; Mtau=params.M_tau; Gnu=params.G_nu;
    Q=Mtau*Gnu; samples=(0:M-1).';
    phases=complex(zeros(M,Mtau,Gnu));
    for delay_index=1:Mtau
        phases(:,delay_index,:)=reshape(exp(1j*2*pi* ...
            (samples-params.delay_grid_2023(delay_index))* ...
            params.doppler_grid_2023(:).'/(M*params.N)),M,1,Gnu);
    end
    xi=complex(zeros(M));
    for transmitter=1:Nt
        transmitter_offset=(transmitter-1)*Q;
        for delay_i=1:Mtau
            delay_i_value=params.delay_grid_2023(delay_i);
            phase_i=reshape(phases(:,delay_i,:),M,Gnu);
            indices_i=transmitter_offset+delay_i:Mtau:transmitter_offset+Q;
            for delay_j=1:Mtau
                delay_j_value=params.delay_grid_2023(delay_j);
                output_j=mod(samples+delay_j_value-delay_i_value,M)+1;
                phase_j=reshape(phases(output_j,delay_j,:),M,Gnu);
                indices_j=transmitter_offset+delay_j:Mtau:transmitter_offset+Q;
                weights=sum((phase_i*covariance(indices_i,indices_j)).*conj(phase_j),2);
                matrix_indices=(1:M).'+(output_j-1)*M;
                xi(matrix_indices)=xi(matrix_indices)+weights;
            end
        end
    end
    xi=(xi+xi')/2;
end
function X=local_affine(symbols,precoder,M,N,Nt)
    X=complex(zeros(M,N,Nt)); for tx=1:Nt, X(:,:,tx)=symbols(:,:,tx)*precoder'; end
end
function [M,N,Nt,Nr,data_columns,mod_order,data_power]=local_validate_inputs(Ydd,Hhat,frame,params,n0,covariance)
    required={'M','N','Nt','Nr','K2','M_tau','G_nu','delay_grid_2023','doppler_grid_2023','ModOrder'};
    if ~isstruct(params)||~all(isfield(params,required)), error('detect_lmmse_2023:MissingParameter','params is missing required detector fields.'); end
    M=params.M;N=params.N;Nt=params.Nt;Nr=params.Nr;data_columns=N-params.K2;mod_order=params.ModOrder;
    if any([M N Nt Nr params.K2 params.M_tau params.G_nu]<1)||params.K2>=N, error('detect_lmmse_2023:InvalidParams','Invalid channel dimensions.'); end
    data_power=1;if isfield(params,'data_power_2023'),data_power=params.data_power_2023;end
    if ~isscalar(data_power)||~isreal(data_power)||~isfinite(data_power)||data_power<=0,error('detect_lmmse_2023:InvalidDataPower','data_power_2023 must be positive.');end
    if ~isnumeric(Ydd)||~isequal(size(Ydd),[M,N,Nr])||any(~isfinite(Ydd(:))),error('detect_lmmse_2023:InvalidObservation','Ydd must be finite M-by-N-by-Nr.');end
    Q=params.M_tau*params.G_nu;if ~isnumeric(Hhat)||~isequal(size(Hhat),[Q*Nt,Nr])||any(~isfinite(Hhat(:))),error('detect_lmmse_2023:InvalidChannel','Hhat has incompatible dimensions.');end
    if ~isstruct(frame)||~isfield(frame,'D')||~isfield(frame,'P')||~isfield(frame,'Xp')||~isfield(frame,'bits')||~isequal(size(frame.D),[N,data_columns])||~isequal(size(frame.P),[N,params.K2])||~isequal(size(frame.Xp),[M,params.K2,Nt]),error('detect_lmmse_2023:InvalidFrame','frame is incompatible.');end
    if ~isscalar(n0)||~isreal(n0)||~isfinite(n0)||n0<0,error('detect_lmmse_2023:InvalidNoiseVariance','n0 must be finite nonnegative.');end
    if ~(isempty(covariance)||isnumeric(covariance)||isstruct(covariance)),error('detect_lmmse_2023:InvalidPosteriorCovariance','posterior covariance must be numeric or factored.');end
end
