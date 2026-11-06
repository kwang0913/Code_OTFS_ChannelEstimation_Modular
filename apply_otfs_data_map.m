function out = apply_otfs_data_map(in, flag, frame)
%APPLY_OTFS_DATA_MAP Linear DD-data-to-padded-time map or its exact adjoint.
%   frame.data_indices selects non-guard DD entries in M-by-N-by-Nt order.
%   frame.reserved_tf masks all TF pilot/guard coordinates on every antenna.
%   Known pilot values are excluded from this linear operator.

    M = frame.M;
    N = frame.N;
    Nt = frame.Nt;
    L = frame.padLen;

    if strcmp(flag, 'notransp')
        dd = zeros(M, N, Nt);
        dd(frame.data_indices) = in;
        tf = otfs_dd_to_tf(dd, M, N);
        tf(frame.reserved_tf) = 0;
        time = ifft(tf, [], 1) * sqrt(M);
        padded = zeros(frame.numSamps, Nt);
        for t = 1:Nt
            padded(:,t) = apply_padding(time(:,:,t), L, frame.padType);
        end
        out = padded(:);
    elseif strcmp(flag, 'transp')
        padded = reshape(in, frame.numSamps, Nt);
        time = zeros(M, N, Nt);
        for t = 1:Nt
            v = padded(:,t);
            switch upper(frame.padType)
                case 'CP'
                    blocks = reshape(v, M + L, N);
                    x = blocks(L+1:end,:);
                    x(M-L+1:M,:) = x(M-L+1:M,:) + blocks(1:L,:);
                case 'RCP'
                    x = v(L+1:end);
                    x(end-L+1:end) = x(end-L+1:end) + v(1:L);
                    x = reshape(x, M, N);
                case 'ZP'
                    blocks = reshape(v, M + L, N);
                    x = blocks(1:M,:);
                case {'RZP', 'NONE'}
                    x = reshape(v(1:M*N), M, N);
                otherwise
                    error('apply_otfs_data_map:UnknownPadding', ...
                        'Unknown padType ''%s''.', frame.padType);
            end
            time(:,:,t) = x;
        end
        tf = fft(time, [], 1) / sqrt(M);
        tf(frame.reserved_tf) = 0;
        dd = otfs_tf_to_dd(tf, M, N);
        out = dd(frame.data_indices);
    else
        error('apply_otfs_data_map:UnknownTranspose', ...
            'Expected ''notransp'' or ''transp'', got ''%s''.', flag);
    end
end
