function [ber,bits_hat,soft,diagnostics] = ...
        detect_lmmse_2022_v2(received,Hhat,frame,params,n0)
%DETECT_LMMSE_2022_V2 Jointly detect the full RCP-protected data frame.
    M = params.M; N = params.N; Q = params.M_tau*params.G_nu;
    sample_count = M*N;
    y = reshape(received(:,frame.data_locations,:),sample_count*params.Nr,1);
    gains = local_gains(Hhat,params,Q);
    forward = @(x) local_apply(x,gains,params);
    adjoint = @(z) local_adjoint(z,gains,params);
    rhs = adjoint(y);
    dimension = sample_count*params.Nt;
    regularization = max(n0/params.data_power_2022,eps);
    if dimension <= 256
        channel = complex(zeros(sample_count*params.Nr,dimension));
        for column = 1:dimension
            basis = complex(zeros(dimension,1));
            basis(column) = 1;
            channel(:,column) = forward(basis);
        end
        normal = channel'*channel+regularization*eye(dimension);
        normal = (normal+normal')/2;
        xhat = normal\rhs;
        iterations = 1;
        relative_residual = norm(normal*xhat-rhs)/max(norm(rhs),eps);
        converged = true;
        preconditioner_name = 'direct';
    else
        inverse_blocks = local_block_jacobi_inverse( ...
            gains,params,regularization);
        preconditioner = @(x) local_apply_block_inverse( ...
            x,inverse_blocks,sample_count,params.Nt);
        preconditioner_name = 'block-jacobi';
        normal_map = @(x) adjoint(forward(x))+regularization*x;
        [xhat,flag,relative_residual,iterations] = pcg(normal_map,rhs, ...
            params.lmmse_tol_2022,params.lmmse_max_iter_2022, ...
            preconditioner);
        converged = flag==0;
        if flag ~= 0
            error('detect_lmmse_2022_v2:NotConverged', ...
                'PCG failed with flag %d and residual %.3e after %d iterations.', ...
                flag,relative_residual,iterations);
        end
    end
    diagnostics = struct('iterations',iterations, ...
        'relative_residual',relative_residual,'converged',converged, ...
        'preconditioner',preconditioner_name);
    time_hat = reshape(xhat,M,N,params.Nt);
    F = fft(eye(N))/sqrt(N);
    soft = complex(zeros(size(time_hat)));
    for tx = 1:params.Nt
        soft(:,:,tx) = time_hat(:,:,tx)*F;
    end
    hard = qamdemod(soft(:)/sqrt(params.data_power_2022), ...
        params.ModOrder,UnitAveragePower=true);
    bits_hat = reshape(int2bit(hard,log2(params.ModOrder)), ...
        size(frame.source_bits));
    ber = mean(bits_hat(:) ~= frame.source_bits(:));
end

function y = local_apply(x,gains,params)
    sample_count = params.M*params.N;
    y = complex(zeros(sample_count,params.Nr));
    for tx = 1:params.Nt
        input = x((tx-1)*sample_count+(1:sample_count));
        for delay_index = 1:params.M_tau
            y = y+gains{delay_index,tx}.*circshift(input, ...
                params.delay_grid_2022(delay_index));
        end
    end
    y = y(:);
end

function x = local_adjoint(y,gains,params)
    sample_count = params.M*params.N;
    y = reshape(y,sample_count,params.Nr);
    x = complex(zeros(sample_count*params.Nt,1));
    for tx = 1:params.Nt
        accumulator = complex(zeros(sample_count,1));
        for delay_index = 1:params.M_tau
            accumulator = accumulator+circshift( ...
                sum(conj(gains{delay_index,tx}).*y,2), ...
                -params.delay_grid_2022(delay_index));
        end
        x((tx-1)*sample_count+(1:sample_count)) = accumulator;
    end
end

function gains = local_gains(Hhat,params,Q)
    phase = local_phase(params);
    gains = cell(params.M_tau,params.Nt);
    for tx = 1:params.Nt
        for delay_index = 1:params.M_tau
            rows = (tx-1)*Q+delay_index:params.M_tau:tx*Q;
            gains{delay_index,tx} = phase*Hhat(rows,:);
        end
    end
end

function inverse_blocks = local_block_jacobi_inverse( ...
        gains,params,regularization)
    sample_count = params.M*params.N;
    Nt = params.Nt;
    blocks = complex(zeros(Nt,Nt,sample_count));
    for tx = 1:Nt
        for other_tx = tx:Nt
            entry = complex(zeros(sample_count,1));
            for delay_index = 1:params.M_tau
                same_output_energy = sum( ...
                    conj(gains{delay_index,tx}).* ...
                    gains{delay_index,other_tx},2);
                entry = entry+circshift(same_output_energy, ...
                    -params.delay_grid_2022(delay_index));
            end
            blocks(tx,other_tx,:) = reshape(entry,1,1,sample_count);
            blocks(other_tx,tx,:) = reshape(conj(entry),1,1,sample_count);
        end
        blocks(tx,tx,:) = blocks(tx,tx,:)+regularization;
    end
    identity_pages = repmat(eye(Nt),1,1,sample_count);
    inverse_blocks = pagemldivide(blocks,identity_pages);
    if any(~isfinite(inverse_blocks(:)))
        error('detect_lmmse_2022_v2:InvalidPreconditioner', ...
            'Block-Jacobi preconditioner contains non-finite values.');
    end
end

function output = local_apply_block_inverse( ...
        input,inverse_blocks,sample_count,Nt)
    input_matrix = reshape(input,sample_count,Nt);
    input_pages = permute(input_matrix,[2 3 1]);
    output_pages = pagemtimes(inverse_blocks,input_pages);
    output_matrix = permute(output_pages,[3 1 2]);
    output = output_matrix(:);
end

function phase = local_phase(params)
    sample_count = params.M*params.N;
    samples = params.pilot_serial_length_2022+params.cp_length_2022+ ...
        (0:sample_count-1).';
    phase = exp(1j*2*pi*samples*params.doppler_grid_2022(:).' / ...
        params.doppler_reference_samples_2022);
end
