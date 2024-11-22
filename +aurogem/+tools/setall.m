function setall(obj, property_suffix, val, opts)
arguments
    obj
    property_suffix (1, :) char
    val
    opts.print (1, 1) logical = false
end
ignored_properties = ["Geoaxes", "Mapaxes"];
property_suffix = char(property_suffix);
l = length(property_suffix);
properties = fieldnames(get(obj, 'factory'));
for n = 1:numel(properties)
    p = properties{n};
    if any(contains(p, ignored_properties))
        continue
    end
    if strcmp(p(end - l + 1:end), property_suffix)
        if opts.print
            fprintf('%s\n', p)
        else
            set(obj, ['default', p(8:end)], val)
        end
    end
end