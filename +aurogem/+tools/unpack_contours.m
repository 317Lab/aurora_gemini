function [x, y] = unpack_contours(cntr,opts)
arguments
    cntr (2,:) double
    opts.rank (1,1) int16 = 1
    opts.min_sep (1,1) double {mustBeNonnegative} = 0
end
%#ok<*AGROW>

mean_ys = [];
lc = size(cntr,2);
i = 1;
j = 1;
while true
    if i >= lc
        break
    end
    N = cntr(2,i);
    x = cntr(1,i+1:i+N);
    y = cntr(2,i+1:i+N);
    i = i+1+N;
    curve = [x; y];
    if abs(curve(1,end) - curve(1,1)) < opts.min_sep
        continue
    end
    mean_ys = [mean_ys,mean(y)];
    curves{j} = curve;
    j = j+1;
end

[~,ids] = sort(mean_ys);
curves_sorted = curves(ids);
if opts.rank == -1
    opts.rank = length(mean_ys);
end
curve = curves_sorted{opts.rank};
x = curve(1,:);
y = curve(2,:);
