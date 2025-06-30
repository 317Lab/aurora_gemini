function setup(data_direc)
arguments
    data_direc (1, :) char {mustBeFile}
end

events = readlines(fullfile(data_direc, 'events.dat'));
events_ids = find(startsWith(events, "------------"));
event_titles = strsplit(events(events_ids(1) + 1));
event = struct;
for e = events(events_ids(1):events_ids(2))'
    tmp = strsplit(e);
    if str2double(tmp{1}) == event_id
        for t = 1:length(event_titles)
            event.(event_titles(t)) = tmp{t};
        end
        break
    end
end

aurogem.glow.generate_invert_table(datetime(2023,2,10,9,51,27), 65.12, 212.8, out_direc='data/swop/dasc/glow')
end