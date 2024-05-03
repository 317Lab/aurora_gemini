function [f1,f2] = nearest_factors(n)
    closest_diff = Inf;
    f1 = 1;
    f2 = n;
    for i = 1:sqrt(double(n))
        if rem(n,i) == 0
            j = n/i;
            diff = abs(i-j);
            if diff < closest_diff
                closest_diff = diff;
                f1 = i;
                f2 = j;
            end
        end
    end
    f1 = double(f1);
    f2 = double(f2);
end