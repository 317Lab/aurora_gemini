function setup(data_direc, opts)
arguments
    data_direc (1, :) char {mustBeFolder}
    opts.glat (1, 1) double = 65.12
    opts.glon (1, 1) double = 212.8
    opts.out_direc (1, :) char = ''
end

if isempty(opts.out_direc)
    out_direc = fullfile(data_direc, 'dasc', 'glow');
else
    out_direc = opts.out_direc;
end
if not(exist('out_direc', 'dir'))
    mkdir(out_direc)
end

events = readlines(fullfile(data_direc, 'events.dat'));
events_ids = find(startsWith(events, "------------"));
for e = events(events_ids(1)+2:events_ids(2)-1)'
    tmp = strsplit(e);
    time = datetime(tmp{2}, 'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss''Z');
    aurogem.glow.generate_invert_table(time, opts.glat, opts.glon, out_direc=out_direc)
end
end