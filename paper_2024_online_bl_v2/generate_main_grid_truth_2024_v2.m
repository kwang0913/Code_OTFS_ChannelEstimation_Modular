function truth = generate_main_grid_truth_2024_v2(params, chanParams)
%GENERATE_MAIN_GRID_TRUTH_2024_V2 Map physical paths to nearest grid coefficients.
% Physical observations retain the unsnapped Doppler values in chanParams.
    required = {'pathGains','pathDelays','pathDopplers','AoAs','AoDs'};
    for index = 1:numel(required)
        if ~isfield(chanParams, required{index})
            error('generate_main_grid_truth_2024_v2:MissingChannelParameter', ...
                'chanParams.%s is required.', required{index});
        end
    end
    count = numel(chanParams.pathDelays);
    Q = params.M_tau*params.G_nu;
    truth.H = complex(zeros(Q*params.Nt, params.Nr));
    truth.support = false(Q, 1);
    truth.pathGains = chanParams.pathGains(:).';
    truth.pathGridDelays = zeros(1, count);
    truth.pathGridDopplers = zeros(1, count);
    truth.pathAtomIndices = zeros(1, count);
    for path = 1:count
        [~, delay_index] = min(abs(params.delay_grid_2024-chanParams.pathDelays(path)));
        [~, doppler_index] = min(abs(params.doppler_grid_2024-chanParams.pathDopplers(path)));
        atom = (doppler_index-1)*params.M_tau+delay_index;
        ar = exp(1j*2*pi*(0:params.Nr-1).' * params.ddlr * ...
            sind(chanParams.AoAs(path)));
        at = exp(1j*2*pi*(0:params.Nt-1).' * params.ddlt * ...
            sind(chanParams.AoDs(path)));
        links = chanParams.pathGains(path)*ar*at';
        for tx = 1:params.Nt
            truth.H((tx-1)*Q+atom, :) = ...
                truth.H((tx-1)*Q+atom, :) + links(:, tx).';
        end
        truth.support(atom) = true;
        truth.pathGridDelays(path) = params.delay_grid_2024(delay_index);
        truth.pathGridDopplers(path) = params.doppler_grid_2024(doppler_index);
        truth.pathAtomIndices(path) = atom;
    end
end
