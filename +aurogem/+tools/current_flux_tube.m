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
    opts.stream3_step (1, 1) double {mustBePositive} = 1e-2
    opts.stream3_maxvert (1, 1) double {mustBePositive} = 2e6
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

% unpack data
jx = double(permute(dat.J2(lbz:ubz, lbx:ubx, lby:uby), [3, 2, 1])) / scl.x^2;
jy = double(permute(dat.J3(lbz:ubz, lbx:ubx, lby:uby), [3, 2, 1])) / scl.x^2;
jz = double(permute(dat.J1(lbz:ubz, lbx:ubx, lby:uby), [3, 2, 1])) / scl.x^2;

% unpack pars
p0 = pars.p0;
r0 = pars.r(1);
r1 = pars.r(2);
v0 = pars.v0;
v1 = pars.v1;
do_reverse = pars.do_reverse;
split_factor = pars.split_factor;
kink_check = pars.kink_check;

perimeter = pi*(3*(r0+r1)-sqrt((3*r0+r1)*(r0+3*r1)));
resolution = perimeter * pars.resolution / 100;

% tolerance values
gap_tolerance = 2 * min([dx; dy; dz]);
split_tube_tolerance = split_factor * pi * r0 * r1 / resolution;
kink_range_deg = wrapTo360(pars.kink_range_deg);

if opts.debug
    fprintf('gap_tolerance = %.3f %s\n', gap_tolerance, unt.x)
    fprintf('split_tube_tolerance = %.3f %s\n', split_tube_tolerance, unt.x)
end

% starting curve
v0 = v0 / norm(v0);
v1 = v1 - dot(v0, v1) * v0;
v1 = v1 / norm(v1);
t = repmat(linspace(0, 1, resolution)', 1, 3);
c0 = p0 + r0*cos(2*pi*t).*v0 + r1*sin(2*pi*t).*v1;

% streamlines
step = opts.stream3_step;
maxvert = opts.stream3_maxvert;
if do_reverse == 2
    fprintf('Running boh directions (experimental).\n')
    verts1 = stream3(Xm, Ym, Zm, -jx, -jy, -jz , c0(:, 1), c0(:, 2), c0(:, 3), [step maxvert]);
    verts2 = stream3(Xm, Ym, Zm,  jx,  jy,  jz , c0(:, 1), c0(:, 2), c0(:, 3), [step maxvert]);
    verts = cellfun(@(v1, v2) [flip(v1); v2], verts1, verts2, 'UniformOutput', 0);
    c0 = cell2mat(cellfun(@(v) v(1, :)', verts, 'UniformOutput', 0))';
else
    dir = 1-2*do_reverse;
    verts = stream3(Xm, Ym, Zm, dir*jx, dir*jy, dir*jz ...
        , c0(:, 1), c0(:, 2), c0(:, 3), [step maxvert]);
end
points = vertcat(verts{:});
c1 = cell2mat(cellfun(@(v) v(end, :)', verts, 'UniformOutput', 0))';

if do_reverse == 1
    c1_tmp = c1;
    c1 = c0;
    c0 = c1_tmp;
end

c0_all = c0;
c1_all = c1;

% split tubes by gaps and kinks
c0 = split_at_gaps(c0, split_tube_tolerance);
c1 = split_at_gaps(c1, split_tube_tolerance);

if kink_check
    c0 = split_at_kinks(c0, kink_range_deg, gap_tolerance, opts.debug);
    c1 = split_at_kinks(c1, kink_range_deg, gap_tolerance, opts.debug);
end

curves = {c0, c1};
fluxes = nan(1, 2);
areas = cell(1, 2);
area_axes = nan(1, 2);

is_split = false(1, 2);
fluxes_split = cell(1, 2);
areas_split = cell(1, 2);
area_axes_split = cell(1, 2);

for i = 1:2 % start and end curves
    c = curves{i};

    if any(isnan(c(:)))
        is_split(i) = true;
        fprintf('Split detected in flux tube at end %i.\n', i)
        [c, shift, nan_ids] = clean_split(c); % partition split curves into cells by consecutive nans
        curves{i} = c;
        verts = circshift(verts, shift); % rotate to match nan_ids
        verts(nan_ids) = {nan(1, 3)}; % remove associated vertices
    end

    % calculate fluxes
    if not(is_split(i))
        [fluxes(i), areas{i}, area_axes(i)] = get_flux(c, jx, jy, jz, dx, dy, dz, Xm, Ym, Zm);
    else
        fluxes_split_tmp = nan(1, length(c));
        areas_split_tmp = cell(1, length(c));
        area_axes_split_tmp = nan(1, length(c));
        for j = 1:length(c)
            [flux_tmp, area_tmp, area_axes_tmp] = get_flux(c{j}, jx, jy, jz, dx, dy, dz, Xm, Ym, Zm);
            fluxes_split_tmp(j) = flux_tmp;
            areas_split_tmp{j} = area_tmp;
            area_axes_split_tmp(j) = area_axes_tmp;
        end
        fluxes_split{i} = fluxes_split_tmp;
        areas_split{i} = areas_split_tmp;
        area_axes_split{i} = area_axes_split_tmp;

        fluxes(i) = sum(fluxes_split_tmp);
        if all(area_axes_split_tmp == area_axes_split_tmp(1))
            area_axes(i) = area_axes_split_tmp(1);
            areas{i} = areas_split_tmp{1};
            for j = 2:length(areas_split_tmp)
                areas{i} = areas{i} | areas_split_tmp{j};
            end
        else
            area_ids = find(area_axes_split_tmp == 3);
            if isempty(area_ids)
                area_axes(i) = area_axes_split_tmp(1);
                areas{i} = areas_split_tmp{1};
            else
                area_axes(i) = 3;
                areas{i} = areas_split_tmp{area_ids(1)};
                for area_id = area_ids(2:end)
                    areas{i} = areas{i} | areas_split_tmp{area_id};
                end
            end
        end
    end
end

% generate tube
tube.vertices = verts;
tube.caps.start = curves{1};
tube.caps.start_all = c0_all;
tube.caps.end = curves{2};
tube.caps.end_all = c1_all;
tube.points = points;
tube.flux.in = fluxes(1);
tube.flux.out = fluxes(2);
tube.flux.in_axis = area_axes(1);
tube.flux.out_axis = area_axes(2);
tube.flux.area.in = areas{1};
tube.flux.area.out = areas{2};
if is_split(1)
    tube.flux.in_split = fluxes_split{1};
    tube.flux.in_split_axes = area_axes_split{1};
    tube.flux.area.in_split = areas_split{1};
end
if is_split(2)
    tube.flux.out_split = fluxes_split{2};
    tube.flux.out_split_axes = area_axes_split{2};
    tube.flux.area.out_split = areas_split{2};
end

end

% functions
function [flux, reg, areas_axis] = get_flux(curve, vx, vy, vz, dx, dy, dz, Xm, Ym, Zm)
flux = 0;
xx = squeeze(Xm(1, :, 1))';
yy = squeeze(Ym(:, 1, 1));
zz = squeeze(Zm(1, 1, :));
reg = false([length(xx), length(yy)]);
areas_axis = find(range(curve) < [min(dx), min(dy), min(dz)] * 4);
if length(areas_axis) == 1
    c_avg = mean(curve(:, areas_axis));
    switch areas_axis
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
    areas_axis = 0;
    warning('Inaccurate influx due to non-orthogonal start/end curve.')
end
end

function c = split_at_gaps(c, tol)
c(abs(vecnorm(diff(c)')) > tol) = nan;
end

function c = split_at_kinks(c, deg, tol, debug)
if all(range(c) > tol) % check for non-flat curve
    fprintf('Possible kink detected.\n')
    nhat = diff(c) ./ vecnorm(diff(c)')';
    c_angles = acosd(dot(nhat(2:end,:), nhat(1:end-1,:), 2));
    kink_ids = find(c_angles > deg(1) & c_angles < deg(2)) + 1;
    lk = length(kink_ids);

    % remove coplanar kinks
    rang = 6;
    test_ids = mod(kink_ids + (-1-rang:rang), length(c)) + 1;
    for i = 1:lk
        if any(range(c(test_ids(i, :), :)) < (tol / 4))
            kink_ids(i) = nan;
        end
    end
    kink_ids(isnan(kink_ids)) = [];

    if debug
        close all
        hold on
        plot3(c(:, 1), c(:, 2), c(:, 3), '-xk')
        plot3(c(1, 1), c(1, 2), c(1, 3), 'og')
        plot3(c(kink_ids, 1), c(kink_ids, 2), c(kink_ids, 3), 'ob', MarkerSize=20)
        for i = 1:lk
            plot3(c(test_ids(i, :), 1), c(test_ids(i, :), 2), c(test_ids(i, :), 3), 'om', MarkerSize=10)
        end
        input('Press any key to continue...')
    end
    c(kink_ids) = nan;
end
end

function [c_out, shift, nan_ids] = clean_split(c_in)
for shift = 0:length(c_in) % rotate until first curve in front
    if not(isnan(c_in(1, 1))) && isnan(c_in(end, 1))
        break
    end
    c_in = circshift(c_in, 1);
end
nan_ids = find(isnan(c_in));
split_ids = [1; nan_ids(diff(nan_ids) > 3 + 1); length(c_in)];
c_out = cell(1, length(split_ids)-1);
j = 1;
for i = 1:length(c_out)
    c_tmp = c_in(split_ids(i):split_ids(i+1), :);
    c_tmp(isnan(c_tmp), :) = [];
    curve_closed = [c_tmp; c_tmp(1, :)];
    if length(c_tmp) > 3
        c_out(j) = mat2cell(curve_closed, length(curve_closed), 3);
        j = j+1;
    else
        c_out(end) = [];
    end
end

end