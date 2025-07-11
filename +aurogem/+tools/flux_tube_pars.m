function [base, pars] = flux_tube_pars(direc)
arguments
    direc (1, :) char {mustBeFolder}
end

fn = fullfile(direc, 'plots3d', 'flux_tube_config.nml');
assert(exist(fn, 'file'), 'Flux tube configuration, %s, file not found.', fn)

base_vars = [ ...
    "limx", "limy", "limz", ...
    "limn", "limj", "zoom", ...
    "panx", "pany", "ar" ...
    ];
base_default = { ...
    nan(1, 2), nan(1, 2), nan(1, 2), ...
    nan(1, 2), nan(1, 2), [1.11, 1, 1.3], ...
    [-0.11, 0.2, 0.2], [0, 0.2, 0.2], nan(1, 3) ...
    };
tube_vars = [ ...
    "p0", "r", "v0", "v1", ...
    "resolution", "color", "do_reverse", "do_projection", ...
    "kink_check", "kink_range_deg", "split_factor", "max_diff_factor" ...
    ];
tube_default = { ...
    nan(1, 3), nan(1, 2), [1, 0, 0], [0, 1, 0], ...
    300, nan(1, 3), 0, 1, ...
    1, [45, 135], 10, 50 ...
    };

for i = 1:length(base_vars)
    base.(base_vars(i)) = base_default{i};
end

lines = readlines(fn)';

base_id = find(strcmp(lines, '&base'));
tube_ids = find(startsWith(lines, '&tube'));
end_ids = find(strcmp(lines, '/'));

i = 1;
while true
    if any(base_id + i == end_ids)
        break
    end
    tmp = strsplit(lines(base_id + i), '=');
    field = strrep(regexprep(tmp(1), '\t', ''), ' ', '');
    data = str2double(strsplit(tmp(2), ','));
    base.(field) = data;
    i = i + 1;
end

for id = tube_ids
    name = char(lines(id));
    name = strrep(regexprep(name(6:end), '\t', ''), ' ', '');
    for i = 1:length(tube_vars)
        pars.(name).(tube_vars(i)) = tube_default{i};
    end
    i = 1;
    while true
        if any(id + i == end_ids)
            break
        end
        tmp = strsplit(lines(id + i), '=');
        field = strrep(regexprep(tmp(1), '\t', ''), ' ', '');
        data = str2double(strsplit(tmp(2), ','));
        pars.(name).(field) = data;
        i = i + 1;
    end
end
