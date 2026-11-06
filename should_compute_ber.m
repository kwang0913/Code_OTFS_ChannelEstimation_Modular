function compute = should_compute_ber(params)
%SHOULD_COMPUTE_BER Decide whether a trial runs data detection.
%   COMPUTE = SHOULD_COMPUTE_BER(PARAMS) returns PARAMS.compute_ber when it
%   is set, otherwise PARAMS.trials == 1: single runs measure BER by default
%   and multi-trial runs skip detection unless compute_ber = true.

    if isfield(params, 'compute_ber') && ~isempty(params.compute_ber)
        compute = logical(params.compute_ber);
        if ~isscalar(compute)
            error('should_compute_ber:InvalidFlag', ...
                'compute_ber must be empty or a logical scalar.');
        end
    else
        compute = params.trials == 1;
    end
end
