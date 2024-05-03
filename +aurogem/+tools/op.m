function outer = op(u,v)
arguments
    u (1,:)
    v (1,:)
end
    lu = length(u);
    lv = length(v);
    outer = nan(lu,lv);
    for i = 1:lu
        for j = 1:lv
            outer(i,j) = u(i)*v(j);
        end
    end
end