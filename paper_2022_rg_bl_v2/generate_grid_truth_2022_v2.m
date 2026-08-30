function truth = generate_grid_truth_2022_v2(params,chanParams)
%GENERATE_GRID_TRUTH_2022_V2 Map physical paths to nearest RG grid atoms.
    Q = params.M_tau*params.G_nu;
    count = numel(chanParams.pathDelays);
    truth.H = complex(zeros(Q*params.Nt,params.Nr));
    truth.support = false(Q,1);
    truth.pathGridDelays = zeros(1,count);
    truth.pathGridDopplers = zeros(1,count);
    truth.pathAtomIndices = zeros(1,count);
    for path = 1:count
        [~,delay_index] = min(abs(params.delay_grid_2022-chanParams.pathDelays(path)));
        [~,doppler_index] = min(abs(params.doppler_grid_2022-chanParams.pathDopplers(path)));
        atom = delay_index+(doppler_index-1)*params.M_tau;
        ar = exp(1j*2*pi*(0:params.Nr-1).'*params.ddlr* ...
            sind(chanParams.AoAs(path)));
        at = exp(1j*2*pi*(0:params.Nt-1).'*params.ddlt* ...
            sind(chanParams.AoDs(path)));
        links = chanParams.pathGains(path)*ar*at';
        for tx = 1:params.Nt
            truth.H((tx-1)*Q+atom,:) = ...
                truth.H((tx-1)*Q+atom,:)+links(:,tx).';
        end
        truth.support(atom) = true;
        truth.pathGridDelays(path) = params.delay_grid_2022(delay_index);
        truth.pathGridDopplers(path) = params.doppler_grid_2022(doppler_index);
        truth.pathAtomIndices(path) = atom;
    end
end
