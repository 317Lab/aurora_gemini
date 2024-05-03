function [peaks,ids] = peak_detect(list,options)
    arguments
        list (1,:) {mustBeNumeric}
        options.num (1,1) int32 {mustBeNonnegative} = 1
        options.threshold (1,1) {mustBeNumeric} = -inf
        options.smoothness (1,1) {mustBeNumeric} = 0
    end

    ids = nan(size(list));
    for i = 2:length(list)-1
        p = list(i-1);
        c = list(i);
        n = list(i+1);
        if c > options.threshold
            if all((c-[p,n]) < options.smoothness*range(list))
                if all([p,n]<=c)
                    ids(i) = i;
                end
            end
        end
    end
    ids(isnan(ids)) = [];
    [peaks,order] = sort(list(ids),'descend');
    ids = ids(order);
    if options.num > 0
        if length(peaks) < options.num
            warning(['Fewer than ',num2str(options.num),' peaks found.'])
        else
            peaks = peaks(1:options.num);
            ids = ids(1:options.num);
        end
    end
end