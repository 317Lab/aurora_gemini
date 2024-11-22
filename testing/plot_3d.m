direc = fullfile('..', 'public_html', 'Gemini3D', 'swop_20230210_35487_A_05_update');

% load fluxtube parameters
[base, tube_pars] = aurogem.tools.flux_tube_pars(direc);

%% plotting parameters
save_plot = 0;
reload_tubes = 1;
reload_grid = 0;
reload_data = 0;
draft = 0;
tube_list = 1:3;
vind = 1;
sffx = ["iso", "side", "top"];
fntn = 'Arial'; % font name
fnts = 11 * 2; % font size
linw = 2; % line width
angl = [[-30, 32]; [-90, 0]; [0, 90]]; % view angle (°)
qntl = 0.95; % colorbar range quantile
pprw = [6.5, 6.5, 6.5] * 2; % paper width (inches)
pprh = [5, 5, 5] * 2; % paper height (inches)
ppro = [[0.6, 0.4, 1, 0.2]; [0, 0, 0, 0]; [0, 0, 0, 0]]; % paper offset, [left, bottom, right, top] (inches)
clbh = 0.43; % colorbar height (relative)
clc0 = [0, 0, 0]; % start curve color (rgb)
clc1 = [0, 0, 1]; % end curve color (rgb)
clef = [0.2, 0.9, 0.6]; % electric field color (rgb)
clfd = [1, 0.4, 1]; % flow data color (rgb)
stlo = 0.3; % flux tube opacity
offs = 0.5; % projection line offset (km)

zoom = base.zoom; % camera zoom
panx = base.panx; % pan right (°)
pany = base.pany; % pan up (°)
sffx = char(sffx(vind));
angl = angl(vind, :);
pprw = pprw(vind);
pprh = pprh(vind);
ppro = ppro(vind, :);
zoom = zoom(vind);
panx = panx(vind);
pany = pany(vind);

%% load simulation structures
sats = strsplit(direc, '_');
sats = sats{5};

cfg = gemini3d.read.config(direc);
if not(exist('xg', 'var')) || reload_grid
    xg = gemini3d.read.grid(direc);
    xg = aurogem.tools.shrink(xg);
end
if not(exist('dat', 'var')) || reload_data
    dat = gemini3d.read.frame(direc, 'time', cfg.times(end) ...
        , 'vars', ["J1", "J2", "J3", "ne", "v2", "v3"]);
end
fprintf('Simulation data loaded.\n')

scl.x = 1e-3;  unt.x = 'km';
scl.e = 1e3;   unt.e = 'mV/m';
scl.f = 1e-3;  unt.f = 'kA';
scl.n = 1e0;   unt.n = 'm^{-3}'; clm.n = 'L9';
scl.j = 1e6;   unt.j = 'uA/m^2'; clm.j = 'D1A';
scl.v = 1e-3;  unt.v = 'km/s';
scl.qv = 1e-2;

lim.x = base.limx;
lim.y = base.limy;
lim.z = base.limz;

% unpack grid
x = double(xg.x2(3:end-2) * scl.x);
y = double(xg.x3(3:end-2) * scl.x);
z = double(xg.x1(3:end-2) * scl.x);
[~, lbx] = min(abs(x - max(min(x), lim.x(1))));
[~, lby] = min(abs(y - max(min(y), lim.y(1))));
[~, lbz] = min(abs(z - max(min(z), lim.z(1))));
[~, ubx] = min(abs(x - min(max(x), lim.x(2))));
[~, uby] = min(abs(y - min(max(y), lim.y(2))));
[~, ubz] = min(abs(z - min(max(z), lim.z(2))));

x = x(lbx:ubx);
y = y(lby:uby);
z = z(lbz:ubz);
dx = xg.dx2h(lbx:ubx) * scl.x;
dy = xg.dx3h(lby:uby) * scl.x;
dz = xg.dx1h(lbz:ubz) * scl.x;
[Xm, Ym, Zm] = meshgrid(x, y, z);
[dXm, dYm] = meshgrid(dx, dy);
Bmag = mean(xg.Bmag(:));

reverse_test_tolerance = min([dx; dy; dz]);
ar = [range(x), range(y), range(z)];

% load track data
track_fn = fullfile(direc, 'ext', 'tracks.h5');
for s = sats
    track_x.(s) = h5read(track_fn, ['/', s, '/Coordinates/Magnetic/East']) * scl.x;
    track_y.(s) = h5read(track_fn, ['/', s, '/Coordinates/Magnetic/North']) * scl.x;
    track_vx.(s) = h5read(track_fn, ['/', s, '/Flow/Magnetic/East'])';
    track_vy.(s) = h5read(track_fn, ['/', s, '/Flow/Magnetic/North'])';
end

% unpack data
j1 = squeeze(dat.J1(ubz, lbx:ubx, lby:uby))'*scl.j; % A/km^2 = uA/m^2
ne = log10(squeeze(dat.ne(lbz:ubz, length(x)/2, lby:uby))')*scl.n; % m^-3
v2 = permute(squeeze(dat.v2(ubz-1:ubz, lbx:ubx, lby:uby)), [3, 2, 1])*scl.v; % km/s
v3 = permute(squeeze(dat.v3(ubz-1:ubz, lbx:ubx, lby:uby)), [3, 2, 1])*scl.v; % km/s
E2 = -v3*Bmag*scl.e/scl.v; % mV/m
E3 = v2*Bmag*scl.e/scl.v; % mV/m
E2(:, :, 1) = nan;
E3(:, :, 1) = nan;

%% generate current flux tubes
if reload_tubes || not(exist('tubes', 'var'))
    for t = fieldnames(tube_pars)'
        tube_name = t{1};
        fprintf('Loading tube %s...\n', tube_name)
        tubes.(tube_name) = aurogem.tools.current_flux_tube(xg, dat, ...
            tube_pars.(tube_name), xlims = lim.x, ylims = lim.y, zlims = lim.z);
    end
end

%% plotting
close all
reset(0)
set(0, 'defaultSurfaceEdgeColor', 'flat')
set(0, 'defaultLineLineWidth', linw)
set(0, 'defaultQuiverLineWidth', linw)
jules.tools.setall(0, 'FontName', fntn)
jules.tools.setall(0, 'FontSize', fnts)
jules.tools.setall(0, 'Multiplier', 1)
colorcet = @jules.tools.colorcet;

lim.j = [-1, 1]*quantile(abs(j1(:)), qntl);
lim.n = [quantile(abs(ne(:)), 1-qntl), quantile(abs(ne(:)), qntl)];

close all
fig = figure;
set(gcf, 'Units', 'inches', 'Position', [0, 0, pprw, pprh])
axj = axes;
axn = axes;
axa = [axj, axn];
set(axj, 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none')
set(axn, 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none')
set(axa, 'units', 'inches', 'position', ...
    [ppro(1), ppro(2), pprw - ppro(1) - ppro(3), pprh - ppro(2) - ppro(4)])
hold on

% current flux tubes
n = 1;
in0 = false(size(j1));
out = false(size(j1));
for t = fieldnames(tube_pars)'
    tube_name = t{1};
    pars = tube_pars.(tube_name);
    color = pars.color;
    tube = tubes.(tube_name);
    verts = tube.vertices;
    points = tube.points;
    c0 = tube.caps.start;
    c1 = tube.caps.end;
    in0 = in0 | tube.flux.area.in;
    out = out | tube.flux.area.out;
    flux0 = tube.flux.in*scl.f;
    flux1 = tube.flux.out*scl.f;
    outline = tube.outline;
    inline = tube.inline;
    fprintf('Tube %s: influx = %.2f %s, outflux = %.2f %s, ratio = %.2f.\n' ...
        , tube_name, flux0, unt.f, flux1, unt.f, flux0 / flux1)
    flux_string_in.(tube_name) = ['\color[rgb]{', num2str(color), '}'...
        , num2str(flux0, '%3.1f'), ' ', unt.f, '  ', '\color{black}'];
    flux_string_out.(tube_name) = ['\color[rgb]{', num2str(color), '}'...
        , num2str(flux1, '%3.1f'), ' ', unt.f, '  ', '\color{black}'];

    stl = streamline(verts);
    plot3(c0(:, 1), c0(:, 2), c0(:, 3), 'Color', clc0);
    if vind == 1
        plot3(c0(:, 1), c0(:, 2), c0(:, 3)*0+z(1)+offs, 'Color', clc0, 'LineStyle', ':');
    end
    if ~iscell(c1)
        plot3(c1(:, 1), c1(:, 2), c1(:, 3), 'Color', clc1);
        if vind == 1
            plot3(c1(:, 1), c1(:, 2), c1(:, 3)*0+z(1)+offs, 'Color', clc1, 'LineStyle', ':');
        end
    else
        for i=1:length(c1)
            c = cell2mat(c1(i));
            plot3(c(:, 1), c(:, 2), c(:, 3), 'Color', clc1);
            if vind == 1
                plot3(c(:, 1), c(:, 2), c(:, 3)*0+z(1), 'Color', clc1, 'LineStyle', ':');
            end
        end
    end
    if vind == 1
        plot3(outline(:, 1)*0+x(end)-offs, outline(:, 2), outline(:, 3), 'Color', color)
    end
    if pars.do_inline && vind == 1
        plot3(inline(:, 1)*0+x(end)-offs, inline(:, 2), inline(:, 3), 'Color', color)
    end
    scatter3(points(:, 1)*0+x(end)-offs, points(:, 2), points(:, 3), 'MarkerEdgeColor', color)
    alpha(0.01)
    set(stl, 'Color', [color, stlo], 'LineWidth', linw/2)
end

% field-aligned current slice
j1(in0 | out) = j1(in0 | out) / 2;
j1_slice = repmat(j1, [1, 1, length(z)]);
slice(axj, Xm, Ym, Zm, -j1_slice, [], [], z(1))
colormap(axj, colorcet(clm.j))
clim(axj, lim.j)
if vind ==1
    clb = colorbar(axj);
    clb.Label.String = sprintf('j_{||} (%s)', unt.j);
    clb.Position = [0.89, clbh+2*(1-2*clbh)/3, 0.015, clbh];
end

% electron density slice
ne_slice = permute(repmat(ne, [1, 1, length(x)]), [1, 3, 2]);
slice(axn, Xm, Ym, Zm, ne_slice, x(end), [], [])
colormap(axn, colorcet(clm.n))
clim(axn, lim.n)
if vind == 1
    clb = colorbar(axn);
    clb.Label.String = sprintf('log_{10} n_e (%s)', unt.n);
    clb.Position = [0.89, (1-2*clbh)/3, 0.015, clbh];
end

% track data
if vind ~= 2
    for s = sats
        x_tmp = track_x.(s);
        y_tmp = track_y.(s);
        z_tmp = ones(size(x_tmp)) * z(1);
        vx_tmp = track_vx.(s) * scl.qv;
        vy_tmp = track_vy.(s) * scl.qv;
        vz_tmp = zeros(size(x_tmp));
        quiver3(x_tmp, y_tmp, z_tmp, vx_tmp, vy_tmp, vz_tmp, 0, '.-', 'Color', clfd)
        quiver3(x(1), y(end), z(end)*0.9, scl.qv, 0, 0, 0, '.-', 'Color', clfd)
    end
end

% electric field
if vind ~= 2
    n_qx = 3;
    n_qy = 4;
    qx_ids = round(linspace(1/2/n_qx, 1-1/2/n_qx, n_qx)*length(x));
    qy_ids = round(linspace(1/2/n_qy, 1-1/2/n_qy, n_qy)*length(y));
    Xm_q = Xm(qy_ids, qx_ids, 1:2);
    Ym_q = Ym(qy_ids, qx_ids, 1:2);
    Zm_q = Zm(qy_ids, qx_ids, 1:2);
    E2_q = E2(qy_ids, qx_ids, :);
    E3_q = E3(qy_ids, qx_ids, :);
    E1_q = zeros(size(E2_q));
    anne = vecnorm([E2_q(1, end, 2), E3_q(1, end, 2), E1_q(1, end, 2)]);
    quiver3(Xm_q, Ym_q, Zm_q, E2_q, E3_q, E1_q, 0.8, 'Color', clef, 'MarkerSize', 2)
    text(Xm_q(1, end, 2)*0.93, Ym_q(1, end, 2)*1.15, Zm_q(1, end, 2)*1.05, ...
        sprintf('%.0f %s', anne, unt.e), 'Color', clef, 'FontSize', fnts*0.7)
end

% annotation
if vind == 1 && length(tube_list) == 3
    annotation('textbox', [0.02, 0.97, 0.01, 0.01], 'String', [...
        'Current In:  ', flux_string_in.(tube_name), newline, ...
        'Current Out:  ', flux_string_out.(tube_name), newline, ...
        ], 'FitBoxToText', 'on', 'EdgeColor', 'none', 'BackgroundColor', 'none', ...
        'FontSize', fnts, 'FontName', fntn)
end

% adjust view positioning
view(axa, angl)
xlim(axa, 1.03*lim.x); ylim(axa, 1.01*lim.y); zlim(axa, lim.z.*[1, 1.01])
pbaspect(axj, ar); pbaspect(axn, ar)

xlabel(axj, sprintf('Mag. east (%s)', unt.x))
ylabel(axj, sprintf('Mag. north (%s)', unt.x))
zlabel(axj, sprintf('Mag. up (%s)', unt.x))
xlabel(axn, sprintf('Mag. east (%s)', unt.x))
ylabel(axn, sprintf('Mag. north (%s)', unt.x))
zlabel(axn, sprintf('Mag. up (%s)', unt.x))

camzoom(axj, zoom); camzoom(axn, zoom)
camdolly(axj, -panx, -pany, 0); camdolly(axn, -panx, -pany, 0)
box(axa, 'on')

% save figure
if save_plot
    [~, filename] = fileparts(direc);
    if ~isempty(sffx)
        filename = [filename, '_', sffx, '.png'];
    else
        filename = [filename, '.png'];
    end
    filename = fullfile(direc, 'plots3d', filename);
    fprintf('Saving %s\n', filename)
    print(fig, filename, '-dpng', '-r96');
    close all
end
