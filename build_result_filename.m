function fname = build_result_filename(params)
%BUILD_RESULT_FILENAME Build a descriptive filename from simulation params.
%   fname = build_result_filename(params)
%
%   Example: 'std_QAM4_pp16_tgt1234_4x16_512x128_pb2d63D63.mat'

    % Config name (std / hf)
    if isfield(params, 'config_name')
        cfg = params.config_name;
    else
        cfg = 'sim';
    end

    % Modulation: QAM4, PSK16, etc.
    if isfield(params, 'modulation')
        mod_type = upper(params.modulation);
    else
        mod_type = 'QAM';
    end
    mod_str = sprintf('%s%d', mod_type, params.ModOrder);

    % Pilot power
    pp_str = sprintf('pp%g', params.pilot_power);

    % Target indices
    tgt_str = sprintf('tgt%s', strrep(num2str(params.selected_targets(:)'), ' ', ''));

    % Antenna config
    ant_str = sprintf('%dx%d', params.Nt, params.Nr);

    % Grid size
    grid_str = sprintf('%dx%d', params.M, params.N);

    % Pilot blocks: e.g. 'pb2d63D63' — lowercase d=delay block, uppercase D=Doppler block, number=pilot count
    if isfield(params, 'pilot_blocks') && ~isempty(params.pilot_blocks)
        pb = params.pilot_blocks;
        parts = cell(1, length(pb));
        for i = 1:length(pb)
            if strcmp(pb(i).dim, 'delay')
                parts{i} = sprintf('d%d', pb(i).n_pilots);
            else
                parts{i} = sprintf('D%d', pb(i).n_pilots);
            end
        end
        pb_str = sprintf('pb%d%s', length(pb), strjoin(parts, ''));
    else
        pb_str = 'pb0';
    end

    ts = datestr(now, 'yyyymmdd_HHMMSS');
    fname = sprintf('%s_%s_%s_%s_%s_%s_%s_%s.mat', cfg, mod_str, pp_str, tgt_str, ant_str, grid_str, pb_str, ts);
end
