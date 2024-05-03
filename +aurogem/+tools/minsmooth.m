function [y,w] = minsmooth(x)
arguments
    x (1,:) {mustBeNumeric}
end
    w = 0;
    is_unique = false;
    if numel(unique(x)) == 1
        y = x;
        return
    end
    if any(diff(x) < 0)
        error('jules:tools:minsmooth:inputError','List not monotonically')
    end
    while not(is_unique)
        w = w+1;
        y = smooth(x,w);
        is_unique = numel(y) == numel(unique(y));
    end
end