function [edges,theta_max] = find_max_edges(array,opts)
    arguments
        array (:,:) double {mustBeNumeric}
        opts.res (1,1) int32 {mustBePositive} = 128
        opts.theta (1,1) double = inf;
    end
    
    Gx = [[-1,0,1];[-2,0,2];[-1,0,1]];
    Gy = [[1,2,1];[0,0,0];[-1,-2,-1]];
    
    conv_x = abs(conv2(Gx,array));
    conv_y = abs(conv2(Gy,array));
    conv_x = conv_x(3:end-2,3:end-2);
    conv_y = conv_y(3:end-2,3:end-2);
    
    if isinf(opts.theta)
        dtheta = 2*pi/double(opts.res);
        thetas = 0:dtheta:pi;
        sums = zeros(size(thetas));
        for itheta = 1:length(thetas)
            theta = thetas(itheta);
            sums(itheta) = sum(cos(theta)*conv_x + sin(theta)*conv_y,'all');
        end
        [~,imax] = max(sums);
        theta_max = thetas(imax);
    else
        theta_max = opts.theta;
    end
    edges = cos(theta_max)*conv_x + sin(theta_max)*conv_y;
    edges = rescale(edges);
end