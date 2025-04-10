function tube = current_flux_tube(xg, dat, pars, opts)
arguments
    xg (1, 1) struct
    dat (1, 1) struct
    pars (1, 1) struct
    opts.do_plot (1, 1) logical = false
    opts.xlims (1, 2) double = [-1, 1]*inf
    opts.ylims (1, 2) double = [-1, 1]*inf
    opts.zlims (1, 2) double = [-1, 1]*inf
    opts.do_reverse_test (1, 1) logical = false
    opts.debug (1, 1) logical = false
    opts.outline_method (1, 1) int32 = 0
end

% unpack grid
scl.x = 1e-3;
unt.x = 'km';
x = double(xg.x2(3:end-2)*scl.x);
y = double(xg.x3(3:end-2)*scl.x);
z = double(xg.x1(3:end-2)*scl.x);
[~, lbx] = min(abs(x-max(min(x), opts.xlims(1))));
[~, lby] = min(abs(y-max(min(y), opts.ylims(1))));
[~, lbz] = min(abs(z-max(min(z), opts.zlims(1))));
[~, ubx] = min(abs(x-min(max(x), opts.xlims(2))));
[~, uby] = min(abs(y-min(max(y), opts.ylims(2))));
[~, ubz] = min(abs(z-min(max(z), opts.zlims(2))));

x = x(lbx:ubx);
y = y(lby:uby);
z = z(lbz:ubz);
dx = xg.dx2h(lbx:ubx)*scl.x;
dy = xg.dx3h(lby:uby)*scl.x;
dz = xg.dx1h(lbz:ubz)*scl.x;
[Xm, Ym, Zm] = meshgrid(x, y, z);
lz = length(z);

% unpack data
jx = double(permute(dat.J2(lbz:ubz, lbx:ubx, lby:uby), [3, 2, 1]))/scl.x^2;
jy = double(permute(dat.J3(lbz:ubz, lbx:ubx, lby:uby), [3, 2, 1]))/scl.x^2;
jz = double(permute(dat.J1(lbz:ubz, lbx:ubx, lby:uby), [3, 2, 1]))/scl.x^2;

% unpack pars
p0 = pars.p0;
r0 = pars.r(1);
r1 = pars.r(2);
v0 = pars.v0;
v1 = pars.v1;
resolution = pars.resolution;
do_reverse = pars.do_reverse;
split_factor = pars.split_factor;
kink_check = pars.kink_check;

% tolerance values
gap_tolerance = 2 * min([dx; dy; dz]);
split_tube_tolerance = split_factor * pi * r0 * r1 / resolution;
max_diff = gap_tolerance * pars.max_diff_factor;
min_kink_deg = 45;

if opts.debug
    fprintf('gap_tolerance = %.3f %s\n', gap_tolerance, unt.x)
    fprintf('split_tube_tolerance = %.3f %s\n', split_tube_tolerance, unt.x)
    fprintf('max_diff = %.3f %s\n', max_diff, unt.x)
end

% starting curve
v0 = v0 / norm(v0);
v1 = v1 - dot(v0, v1) * v0;
v1 = v1 / norm(v1);
t = repmat(linspace(0, 1, resolution)', 1, 3);
c0 = p0 + r0*cos(2*pi*t).*v0 + r1*sin(2*pi*t).*v1;

% streamlines
step = 1e-2;
maxvert = 2e6;

if do_reverse == 2
    fprintf('Running boh directions (experimental).\n')
    verts1 = stream3(Xm, Ym, Zm, jx, jy, jz , c0(:, 1), c0(:, 2), c0(:, 3), [step maxvert]);
    verts2 = stream3(Xm, Ym, Zm, -jx, -jy, -jz , c0(:, 1), c0(:, 2), c0(:, 3), [step maxvert]);
    verts = cellfun(@(v1, v2) [flip(v1); v2], verts1, verts2, 'UniformOutput', 0);
    c0 = cell2mat(cellfun(@(v) v(1, :)', verts, 'UniformOutput', 0))';
else
    dir = 1-2*do_reverse;
    verts = stream3(Xm, Ym, Zm, dir*jx, dir*jy, dir*jz ...
        , c0(:, 1), c0(:, 2), c0(:, 3), [step maxvert]);
end
points = vertcat(verts{:});
c1 = cell2mat(cellfun(@(v) v(end, :)', verts, 'UniformOutput', 0))';
c1_all = c1;

% reverse test
if opts.do_reverse_test && do_reverse ~=2
    for cid = 1:length(c0)
        maxvert = length(cell2mat(verts(cid)));
        vert = cell2mat(stream3(Xm, Ym, Zm, -dir*jx, -dir*jy, -dir*jz ...
            , c1(cid, 1), c1 ...
            (cid, 2), c1(cid, 3), [step maxvert+4]));
        delta = norm(abs(vert(end, :) - c0(cid, :)));
        if delta > gap_tolerance
            warning('Reverse test failed.')
            break
        end
    end
end

% check for split tube
c1(abs(vecnorm(diff(c1)')) > split_tube_tolerance) = nan;

if kink_check && all(range(c1) > gap_tolerance) % check for non-flat curve
    fprintf('Kink detected.\n')
    c1_nhat = diff(c1) ./ vecnorm(diff(c1)')';
    c1_angles = acosd(dot(c1_nhat(2:end,:), c1_nhat(1:end-1,:), 2));
    c1(circshift(c1_angles > min_kink_deg, 1)) = nan; % check for sharp angles
end

if ~any(isnan(c1))
    split = false;
else
    split = true;
    fprintf('Split flux tube detected.\n')
    for shift = 0:length(c1) % rotate until first curve in front
        if not(isnan(c1(1, 1))) && isnan(c1(end, 1))
            break
        end
        c1 = circshift(c1, 1);
    end
    nan_ids = find(isnan(c1));
    split_ids = [1; nan_ids(diff(nan_ids)>3+1); resolution];
    c1_tmp = cell(1, length(split_ids)-1);
    j = 1;
    for i = 1:length(c1_tmp)
        c = c1(split_ids(i):split_ids(i+1), :);
        c(isnan(c), :) = [];
        curve_closed = [c; c(1, :)];
        if length(c) > 10
            c1_tmp(j) = mat2cell(curve_closed, length(curve_closed), 3);
            j = j+1;
        else
            c1_tmp(end) = [];
        end
    end
    c1 = c1_tmp;
    verts = circshift(verts, shift);
    verts(nan_ids) = [];
end

% calculate fluxes
[flux0, in, fid0] = get_flux(c0, jx, jy, jz, dx, dy, dz, Xm, Ym, Zm);
if split
    flux1_split = nan(1, length(c1));
    fid1_split = nan(1, length(c1));
    out_split = cell(1, length(c1));
    for i = 1:length(c1)
        c = cell2mat(c1(i));
        [flux1_tmp, out_tmp, fid1_tmp] = get_flux(c, jx, jy, jz, dx, dy, dz, Xm, Ym, Zm);
        flux1_split(i) = flux1_tmp;
        fid1_split(i) = fid1_tmp;
        out_split{i} = out_tmp;
    end
    flux1 = sum(flux1_split);
    if all(fid1_split == fid1_split(1))
        fid1 = fid1_split(1);
        out = out_split{1};
        for j = 2:length(out_split)
            out = out | out_split{i};
        end
    else
        fid1 = fid1_split(1);
        out = out_split{1};
    end
else
    [flux1, out, fid1] = get_flux(c1, jx, jy, jz, dx, dy, dz, Xm, Ym, Zm);
end

% generate outline and inline
if opts.outline_method == 0
    outline = nan(1, 3);
    inline = nan(1, 3);
elseif opts.outline_method == 1
    ax_id = pars.outline_axis;
    if ax_id == 3
        ax = z;
    elseif ax_id == 2
        ax = y;
    else
        error('The parameter "outline_axis" must be either 2 or 3.')
    end
    ax = ax(1:pars.outline_res:end);
    lax = length(ax);
    outline = nan(2*lax-4, 3);
    inline = nan(2*lax-4, 3);
    for id = 1:lax-2
        ax0 = ax(id);
        ax1 = ax(id+1);
        layer = points(points(:, ax_id) > ax0 & points(:, ax_id) < ax1, :);
        if size(layer, 1) ~= 0
            [~, min_id] = min(layer(:, 5-ax_id));
            [~, max_id] = max(layer(:, 5-ax_id));
            outline(lax-id-1, :) = layer(max_id, :);
            outline(lax-2+id, :) = layer(min_id, :);

            layer_sorted = sort(layer', 5-ax_id)'; %#ok<UDIM>
            [gap, in_id] = max(diff(layer_sorted(:, 5-ax_id)));
            if gap > gap_tolerance
                inline(lax-id-1, :) = layer_sorted(in_id, :);
                inline(lax-2+id, :) = layer_sorted(in_id+1, :);
            end
        end
    end
    outline = outline(not(isnan(outline(:, 1))), :);
    inline = inline(not(isnan(inline(:, 1))), :);

    % clean out/inlines
    break_ids = find(vecnorm(diff(inline)') > max_diff);
    if length(break_ids) > 1
        bad = break_ids(1):break_ids(2);
        inline(bad) = nan;
    end
    break_ids = find(vecnorm(diff(outline)') > max_diff);
    if length(break_ids) > 1
        bad = break_ids(1):break_ids(2);
        outline(bad) = nan;
    end
    outline(vecnorm(diff(outline(:, 2:3))') > max_diff) = nan;

    outline = outline(not(isnan(outline(:, 1))), :);
    inline = inline(not(isnan(inline(:, 1))), :);
elseif opts.outline_method == 2
    [~, id0] = max(cell2mat(cellfun(@(v) mean(v(:, 3)), verts, 'UniformOutput', 0)));
    [~, id1] = min(cell2mat(cellfun(@(v) mean(v(:, 3)), verts, 'UniformOutput', 0)));
    inline = verts{id0};
    outline = verts{id1};
elseif opts.outline_method == 3
    gc = mean(points(:, 2:3));
    vns = cell2mat(cellfun(@(v) mean(vecnorm(v(:, 2:3)' - gc')), verts, 'UniformOutput', 0));
    [~, id0] = min(vns);
    [~, id1] = max(vns);
    inline = verts{id0};
    outline = verts{id1};
else
    error('Outline method not known.')
end

% generate tube
tube.vertices = verts;
tube.caps.start = c0;
tube.caps.end = c1;
tube.caps.end_all = c1_all;
tube.points = points;
tube.flux.in = flux0;
tube.flux.out = flux1;
tube.flux.in_axis = fid0;
tube.flux.out_axis = fid1;
tube.flux.area.in = in;
tube.flux.area.out = out;
tube.outline = outline;
tube.inline = inline;
if split
    tube.flux.out_split = flux1_split;
    tube.flux.out_axis_split = fid1_split;
    tube.flux.area.out_split = out_split;
end

% plotting
if opts.do_plot
    close all
    reset(0)
    set(0, 'defaultLineLineWidth', 2)

    fac = -squeeze(jz(:, :, end));
    fac(in | out) = nan;
    fac = repmat(fac, [1, 1, lz]);

    colorcet = @aurogem.tools.colorcet;
    lim.x = [min(x), max(x)];
    lim.y = [min(y), max(y)];
    lim.z = [min(z), max(z)];
    lim.j = [-1, 1]*quantile(abs(fac(:)), 0.95);

    hold on
    title(sprintf('flux in/out ratio = %f', flux0/flux1))
    stl = streamline(verts);
    plot3(c0(:, 1), c0(:, 2), c0(:, 3), 'g')
    plot3(c0(:, 1), c0(:, 2), ones(size(c0(:, 3)))*lim.z(1), ':g')
    if ~split
        plot3(c1(:, 1), c1(:, 2), c1(:, 3), 'r')
        plot3(c1(:, 1), c1(:, 2), ones(size(c1(:, 3)))*lim.z(1), ':r')
    else
        for i=1:length(c1)
            c = cell2mat(c1(i));
            plot3(c(:, 1), c(:, 2), c(:, 3), 'r')
            plot3(c(:, 1), c(:, 2), ones(size(c(:, 3)))*lim.z(1), ':r')
        end
    end
    slice(Xm, Ym, Zm, fac, [], [], lim.z(1))
    shading flat
    colormap(colorcet('D1A'))
    xlim(lim.x); ylim(lim.y); zlim(lim.z); clim(lim.j)
    view([45, 30])
    pbaspect([range(x), range(y), range(z)])
    set(stl, 'Color', [0, 0, 0, 0.1], 'LineWidth', 0.5)
end

% calculate flux function
    function [flux, reg, fid] = get_flux(curve, vx, vy, vz, dx, dy, dz, Xm, Ym, Zm)
        flux = 0;
        xx = squeeze(Xm(1, :, 1))';
        yy = squeeze(Ym(:, 1, 1));
        zz = squeeze(Zm(1, 1, :));
        reg = false([length(xx), length(yy)]);
        fid = find(range(curve) < [min(dx), min(dy), min(dz)]);
        if length(fid) == 1
            c_avg = mean(curve(:, fid));
            switch fid
                case 1
                    [dYm, dZm] = meshgrid(dy, dz);
                    dAx = (dYm.*dZm)';
                    [~, c_id] = min(abs(c_avg - xx));
                    reg = squeeze(inpolygon(Ym(:, 1, :), Zm(:, 1, :), curve(:, 2), curve(:, 3)));
                    flux = abs(sum(squeeze(vx(:, c_id, :)).*dAx.*reg, 'all'));
                case 2
                    [dXm, dZm] = meshgrid(dx, dz);
                    dAy = (dXm.*dZm)';
                    [~, c_id] = min(abs(c_avg - yy));
                    reg = squeeze(inpolygon(Xm(1, :, :), Zm(1, :, :), curve(:, 1), curve(:, 3)));
                    flux = abs(sum(squeeze(vy(c_id, :, :)).*dAy.*reg, 'all'));
                case 3
                    [dXm, dYm] = meshgrid(dx, dy);
                    dAz = dXm.*dYm;
                    [~, c_id] = min(abs(c_avg - zz));
                    reg = squeeze(inpolygon(Xm(:, :, 1), Ym(:, :, 1), curve(:, 1), curve(:, 2)));
                    flux = abs(sum(squeeze(vz(:, :, c_id)).*dAz.*reg, 'all'));
            end
        else
            fid = 0;
            warning('Inaccurate influx due to non-orthogonal start/end curve.')
        end
    end
end