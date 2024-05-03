function find_factory_properties(s)
arguments
    s (1,:) char {mustBeNonempty}
end
properties = fieldnames(get(0,'factory'));
for n = 1:length(properties)
    p = properties{n};
    if contains(p,s)
            disp(['default',p(8:end)])
    end
end
end
