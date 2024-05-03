function [x2_bounds_A,x3_bounds_A,x2_bounds_B,x3_bounds_B] = find_boundaries(arc,x2_imag,x3_imag,opts)
arguments
    arc
    x2_imag
    x3_imag
    opts.boundary_smoothing_window
    opts.edge_method
    opts.edge_id
end
[X2_imag,X3_imag] = ndgrid(x2_imag,x3_imag);

num_bounds = 2;
edges = tools.find_max_edges(arc,theta=0);
bsw = opts.boundary_smoothing_window;
% fprintf('Boundary smoothing window is approximatly %.0f meters.\n', ...
%     mean(dx3)*double(bsw))
if strcmp(opts.edge_method,'sobel')
    bounds = nan(num_bounds,size(edges,1));
    for i = 1:size(edges,1)
        [~,bounds(:,i)] = tools.peak_detect(edges(i,:), ...
            num=num_bounds,smoothness=0.009);
    end
    x2_bounds_A = sort(x2_imag(2:end-1));
    x3_bounds_A = smoothdata(x3_imag(bounds(1,:)+1),2,'gaussian',bsw);
    x2_bounds_B = x2_bounds_A;
    x3_bounds_B = smoothdata(x3_imag(bounds(2,:)+1),2,'gaussian',bsw);
elseif strcmp(opts.edge_method,'contour')
    edge_id = opts.edge_id; %%% FIX THIS, NEEDS TO BE BETTER SOMEHOW
    [~,cntr_ids] = tools.peak_detect(edges(edge_id,:), ...
        num=num_bounds,smoothness=0.009);
    cntr_vals = arc(edge_id,cntr_ids);
    cntr_A = contour(X2_imag,X3_imag,arc,[1,1]*cntr_vals(1));
    cntr_B = contour(X2_imag,X3_imag,arc,[1,1]*cntr_vals(2));
    close(gcf)
    lc_A = cntr_A(2,1);
    lc_B = cntr_B(2,1);
    fprintf('Primary contour line at %.2f %s\n', ...
        cntr_A(1,1)^ap*scl.arc,unt.arc)
    fprintf('Secondary contour line at %.2f %s\n', ...
        cntr_B(1,1)^ap*scl.arc,unt.arc)
    x2_bounds_A = smoothdata(cntr_A(1,lc_A+3:end),"gaussian");
    x3_bounds_A = smoothdata(cntr_A(2,lc_A+3:end),"gaussian",bsw);
    x2_bounds_B = smoothdata(cntr_B(1,2:lc_B+1),"gaussian");
    x3_bounds_B = smoothdata(cntr_B(2,2:lc_B+1),"gaussian",bsw);
    [~,sort_ids_A] = sort(x2_bounds_A);
    [~,sort_ids_B] = sort(x2_bounds_B);
    x2_bounds_A = x2_bounds_A(sort_ids_A);
    x3_bounds_A = x3_bounds_A(sort_ids_A);
    x2_bounds_B = x2_bounds_B(sort_ids_B);
    x3_bounds_B = x3_bounds_B(sort_ids_B);
end
end