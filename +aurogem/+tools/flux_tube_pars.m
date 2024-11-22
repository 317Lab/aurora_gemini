function [base, pars] = flux_tube_pars(direc)
arguments
    direc (1, :) char {mustBeFolder}
end

fn = fullfile(direc, 'plots3d', 'flux_tube_config.nml');
assert(exist(fn, 'file'), 'Flux tube configuration, %s, file not found.', fn)

base_vars = ["limx", "limy", "limz", "zoom", "panx", "pany"];
base_default = {nan(1, 2), nan(1, 2), nan(1, 2), ones(1, 3), zeros(1, 3), zeros(1, 3)};
tube_vars = ["p0", "r", "v0", "v1", "resolution", "color", "do_reverse", ...
    "do_inline", "split_factor", "max_diff_factor", "outline_axis", "outline_res"];
tube_default = {nan(1, 3), nan(1, 2), [1, 0, 0], [0, 1, 0], 300, nan(1, 3), 0, 1, 10, 100, 3, 2};

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
