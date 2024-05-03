function f = common_factors(lx2,lx3)
    f2 = factor(lx2);
    f3 = factor(lx3);
    p2 = perms(f2)';
    p3 = perms(f3)';
    q2 = unique(p2(1,:));
    q3 = unique(p3(1,:));
    for i = 2:length(f2)-1
        q2 = [q2,unique(prod(p2(1:i,:),1))]; %#ok<AGROW> 
    end
    for i = 2:length(f3)-1
        q3 = [q3,unique(prod(p3(1:i,:),1))]; %#ok<AGROW> 
    end
    f = intersect(q2,q3);
end
