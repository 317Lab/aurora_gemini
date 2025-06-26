function swarm_cdf_to_h5(data_direc)
arguments
    data_direc (1, :) char {mustBeFolder}
end

cdfs = {dir(fullfile(data_direc, '*.cdf')).name};
assert(~isempty(cdfs), 'No cdfs found in %s.', data_direc)
for cdf = cdfs
    is_tii = false;
    is_lps = false;
    if contains(cdf, 'EXPT_EFI')
        is_tii = true;
    elseif contains(cdf, 'OPER_EFI')
        is_lps = true;
    elseif not(contains(cdf, 'OPER_FAC'))
        warning('Skipped %s', cdf{1})
        continue
    end
    path_cdf = fullfile(data_direc, cdf{1});
    path_h5 = [path_cdf(1:end-3), 'h5'];
    [all_data, info] = cdfread(path_cdf, 'DatetimeType', 'datetime');
    vars = info.Variables(:, 1);

    units = info.VariableAttributes.UNITS(:, 2);
    if is_tii
        types = info.VariableAttributes.Type(:, 2);
        descs = info.VariableAttributes.CATDESC(:, 2);
    else
        descs = info.VariableAttributes.DESCRIPTION(:, 2);
    end

    for i = 1:length(vars)
        data = [all_data{:, i}];
        var = ['/', vars{i}];
        unit = units{i};
        desc = descs{i};
        if is_tii
            type = types{i};
            type = lower(type(5:end));
        elseif isdatetime(data)
            type = 'epoch';
        else
            type = 'float';
        end
        if strcmp(type, 'float')
            type = 'double';
        elseif strcmp(type(2:end-1), "int")
            type = [type(1:end-1), '32'];
        elseif strcmp(type, 'epoch')
            type = 'double';
            data = posixtime(data);
            unit = 's';
            desc = 'Unix time.';
        end
        try
            h5create(path_h5, var, size(data), 'Datatype', type)
            h5write(path_h5, var, data)
            h5writeatt(path_h5, var, 'Units', unit)
            h5writeatt(path_h5, var, 'Description', desc)
        catch
            warning('%s in %s already exists', var(2:end), path_h5)
        end

        if strcmp(var, '/Latitude'); lat_gc = double(data); end
        if strcmp(var, '/Radius'); rad = double(data); end
    end

    global_atts = fieldnames(info.GlobalAttributes)';
    for att = global_atts
        h5writeatt(path_h5, '/', att{1}, info.GlobalAttributes.(att{1}))
    end

    if not(is_lps)
        [lat_gd, alt_gd] = geoc2geod(lat_gc, rad);
        try
            h5create(path_h5, '/GeodeticLatitude', size(lat_gd), 'Datatype', 'double')
            h5write(path_h5, '/GeodeticLatitude', lat_gd)
            h5writeatt(path_h5, '/GeodeticLatitude', 'Units', 'degrees')
            h5writeatt(path_h5, '/GeodeticLatitude', 'Description', 'Geodetic Latitude.')
            h5create(path_h5, '/GeodeticAltitude', size(alt_gd), 'Datatype', 'double')
            h5write(path_h5, '/GeodeticAltitude', alt_gd)
            h5writeatt(path_h5, '/GeodeticAltitude', 'Units', 'm')
            h5writeatt(path_h5, '/GeodeticAltitude', 'Description', 'Geodetic Altitude.')
        catch
            warning('Geodetic data already exists.')
        end
    end
end
end