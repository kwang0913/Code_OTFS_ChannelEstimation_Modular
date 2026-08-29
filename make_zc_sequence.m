function seq = make_zc_sequence(n_seq, root_idx)
%MAKE_ZC_SEQUENCE Generate a unit-amplitude Zadoff-Chu sequence.
%   seq = MAKE_ZC_SEQUENCE(n_seq, root_idx) returns a 1 x n_seq CAZAC
%   sequence. root_idx must be coprime with n_seq.

    if n_seq <= 0 || floor(n_seq) ~= n_seq
        error('make_zc_sequence:InvalidLength', 'Sequence length must be a positive integer.');
    end
    if nargin < 2 || isempty(root_idx)
        root_idx = 1;
    end
    root_idx = mod(root_idx, n_seq);
    if root_idx == 0
        root_idx = 1;
    end
    if gcd(root_idx, n_seq) ~= 1
        error('make_zc_sequence:InvalidRoot', 'ZC root must be coprime with sequence length.');
    end

    q = 0:(n_seq - 1);
    if mod(n_seq, 2) == 0
        phase_arg = q.^2;
    else
        phase_arg = q .* (q + 1);
    end
    seq = exp(-1j * pi * root_idx * phase_arg / n_seq);
end
