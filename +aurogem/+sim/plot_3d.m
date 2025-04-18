% Description:
%   Make 3D visualizations of GEMINI simulations including top-boundary FAC
%   slice, central density slice, sparse electric field vectors, track data, and
%   current flux tubes. Tubes parameters are configured in
%   <direc>/plots3d/flux_tube_config.nml (see README.md for details).
%
% Contact:
%   jules.van.irsel.gr@dartmouth.edu
%
% Revisions:
%   02/28/2025  initial implementation (jvi)
%

%#ok<*UNRCH>

% simulations directories
direc_root = fullfile('..', 'public_html', 'Gemini3D');
% direc = fullfile(direc_root, 'swop_20230210_35487_AC_09_SD');
direc = fullfile(direc_root, 'swop_20230210_35487_AC_09_SD');
% direc = fullfile(direc_root, 'swop_20230210_35487_AC_09_unacc_SD');
% direc = fullfile(direc_root, 'swop_20230210_35487_A_09_nobg');
direc_compare = ''; % when comparing tubes of two simulations

% plotting parameters
vinds = 1:3; % view angles to plot (1=iso, 2=side, 3=top)
spins = 0; %0:5:360; % for spinning animation (°)
save_plot = true;
reload_tubes = true; % when only changing plotting parameters
reload_grid = false; % when only changing simulation data
reload_data = false; % when only changing tube configuration
publish_format = 'paper'; % paper or poster
tube_list = 1:3;
debug = 0;

for vind = vinds
for spin = spins
angl = [[-30, 32]; [-90, 0]; [0, 90]]; % view angle (°)
sffx = ["iso", "side", "top"];
fntn = 'Arial'; % font name
qntl = 0.95; % colorbar range quantile
clbh = 0.43; % colorbar height (relative)
clc0 = [0.0, 0.0, 0.0]; % start curve color (rgb)
clc1 = [0.0, 0.0, 0.8]; % end curve color (rgb)
clef = [0.0, 1.0, 1.0]; % electric field color (rgb)
cleb = [1.0, 1.0, 0.0]; % electric field background (rgb)
cltd = [1.0, 0.4, 1.0]; % track data color (rgb)
stlo = 0.3; % flux tube opacity
offs = 0.5; % projection line offset (km)
idef = 3; % electric field legend vind (1, 3)
lnef = 20; % length of electric field background vector (km)
track_data_type = 'current'; % type of track data to plot

if strcmp(publish_format, 'paper')
    fnts = 10 * 2; % font size
    linw = 1.5; % line width
    pprw = [6.5, 2.58, 3.92] * 2; % paper width (inches)
    pprh = [5, 3.02, 3.02] * 2; % paper height (inches)
    clbx = 0.9; % colorbar horizontal position
    clbg = [255, 255, 255] / 255; % background color (rgb)
    % clbg = [221, 221, 221] / 255;
    % clbg = [54, 54, 54] / 255;
    cltx = [0, 0, 0]; % text color (rgb)
    % cltx = [1, 1, 1];
    xrot = [18, 0, 0]; % x label rotation (deg)
    yrot = [-45, 0, 90]; % y label rotation (deg)
    sffx = sffx + "-p";
else
    fnts = 18 * 2; % font size
    linw = 2; % line width
    pprw = [8.5, 4, 4] * 2; % paper width (inches)
    pprh = [7, 4, 3] * 2; % paper height (inches)
    clbx = 0.87; % colorbar horizontal position
    clbg = [20, 21, 20] / 255; % background color (rgb)
    cltx = [1, 1, 1]; % text color (rgb)
    xrot = [18, 0, 0]; % x label rotation (deg)
    yrot = [-45, 0, 90]; % y label rotation (deg)
    % sffx = sffx + "-test";
end

% compare tubes
if exist(fullfile(direc_compare, 'config.nml'), 'file')
    do_compare = true;
else
    do_compare = false;
end

% load fluxtube parameters
if reload_tubes || not(exist('tube_pars', 'var'))
    [base, tube_pars] = aurogem.tools.flux_tube_pars(direc);
end

% animation sping angle
angl(:, 1) = angl(:, 1) + spin;

% view angles
zoom = base.zoom; % camera zoom
panx = base.panx; % pan right (°)
pany = base.pany; % pan up (°)
sffx_tmp = char(sffx(vind));
angl = angl(vind, :);
pprw = pprw(vind);
pprh = pprh(vind);
zoom = zoom(vind);
panx = panx(vind);
pany = pany(vind);
xrot = xrot(vind);
yrot = yrot(vind);

if do_compare
    sffx = [sffx_tmp, '_compare'];
elseif exist('spin', 'var')
    sffx = [sffx_tmp, '_', num2str(spin)];
else
    sffx = sffx_tmp;
end

tube_names = fieldnames(tube_pars);
ntubes = length(fieldnames(tube_pars));
tube_pars_tmp = struct;
for t = tube_list
    tube_pars_tmp.(tube_names{t}) = tube_pars.(tube_names{t});
end
tube_pars = tube_pars_tmp;

%% load simulation structures
sats = strsplit(direc, '_');
sats = sats{5};

cfg = gemini3d.read.config(direc);
if not(exist('xg', 'var')) || reload_grid
    xg = gemini3d.read.grid(direc);
    xg = aurogem.tools.shrink(xg);
end
fprintf('Simulation grid loaded.\n')

if not(exist('dat', 'var')) || reload_data
    dat = gemini3d.read.frame(direc, 'time', cfg.times(end) ...
        , 'vars', ["J1", "J2", "J3", "ne", "v2", "v3"]);
end
if do_compare && (not(exist('dat_compare', 'var')) || reload_data)
    dat_compare = gemini3d.read.frame(direc_compare, 'time', cfg.times(end) ...
        , 'vars', ["J1", "J2", "J3", "ne", "v2", "v3"]);
end
fprintf('Simulation data loaded.\n')

scl.x = 1e-3;  unt.x = 'km';
scl.e = 1e3;   unt.e = 'mV/m';
scl.qe = 2;
scl.f = 1e-3;  unt.f = 'kA';
scl.n = 1e0;   unt.n = 'm^{-3}'; clm.n = 'L9';
scl.j = 1e6;   unt.j = 'uA/m^2'; clm.j = 'D1A';
scl.v = 1e-3;  unt.v = 'km/s';
scl.qv = 1e-2;

lim.x = base.limx;
lim.y = base.limy;
lim.z = base.limz;
lim.n = base.limn;
lim.j = base.limj;

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

if any(isnan(base.ar))
    ar = [range(x), range(y), range(z)];
else
    ar = base.ar;
end

% load track data
track_fn = fullfile(direc, 'ext', 'tracks.h5');
for s = sats
    track_x.(s) = h5read(track_fn, ['/', s, '/Coordinates/Magnetic/East']) * scl.x;
    track_y.(s) = h5read(track_fn, ['/', s, '/Coordinates/Magnetic/North']) * scl.x;
    track_vx.(s) = h5read(track_fn, ['/', s, '/Flow/Magnetic/East'])';
    track_vy.(s) = h5read(track_fn, ['/', s, '/Flow/Magnetic/North'])';
    track_fac.(s) = h5read(track_fn, ['/', s, '/Current/FieldAligned'])';
end

% unpack data
j1 = squeeze(dat.J1(ubz, lbx:ubx, lby:uby))'*scl.j; % A/km^2 = uA/m^2
j2 = squeeze(dat.J2(lbz:ubz, length(x)/2, lby:uby))'*scl.j; % A/km^2 = uA/m^2
ne = log10(squeeze(dat.ne(lbz:ubz, length(x)/2, lby:uby))')*scl.n; % m^-3
v2 = permute(squeeze(dat.v2(ubz-1:ubz, lbx:ubx, lby:uby)), [3, 2, 1])*scl.v; % km/s
v3 = permute(squeeze(dat.v3(ubz-1:ubz, lbx:ubx, lby:uby)), [3, 2, 1])*scl.v; % km/s
E2 = -v3*Bmag*scl.e/scl.v; % mV/m
E3 = v2*Bmag*scl.e/scl.v; % mV/m
E2(:, :, 1) = nan;
E3(:, :, 1) = nan;

% unpack background electric field data
E_bg_filename = fullfile(direc, cfg.E0_dir, gemini3d.datelab(cfg.times(end)) + '.h5');
E2_bg = h5read(E_bg_filename, '/Exit');
E3_bg = h5read(E_bg_filename, '/Eyit');
E2_bg = median(E2_bg(:)) * scl.e;
E3_bg = median(E3_bg(:)) * scl.e;
scl.qe = lnef / vecnorm([E2_bg, E3_bg]);

% generate current flux tubes
if reload_tubes || not(exist('tubes', 'var'))
    for tp = fieldnames(tube_pars)'
        tube_name = tp{1};
        fprintf('Loading tube %s...\n', tube_name)
        tubes.(tube_name) = aurogem.tools.current_flux_tube(xg, dat, ...
            tube_pars.(tube_name), xlims = lim.x, ylims = lim.y, zlims = lim.z, debug=debug);
        if do_compare
            tubes_compare.(tube_name) = aurogem.tools.current_flux_tube(xg, dat_compare, ...
                tube_pars.(tube_name), xlims = lim.x, ylims = lim.y, zlims = lim.z);
        end
    end
end

%% plotting
close all
reset(0)
set(0, 'defaultSurfaceEdgeColor', 'flat')
set(0, 'defaultLineLineWidth', linw)
set(0, 'defaultQuiverLineWidth', linw)
aurogem.tools.setall(0, 'FontName', fntn)
aurogem.tools.setall(0, 'FontSize', fnts)
aurogem.tools.setall(0, 'Multiplier', 1)
colorcet = @aurogem.tools.colorcet;

if any(isnan(lim.j))
    lim.j = [-1, 1]*quantile(abs(j1(:)), qntl);
end
if any(isnan(lim.n))
    lim.n = [quantile(abs(ne(:)), 1-qntl), quantile(abs(ne(:)), qntl)];
end

close all
fig = figure;
set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, pprw, pprh], ...
    'Position', [10, 30, 900, 900*pprh/pprw])
axj = axes;
axn = axes;
axa = [axj, axn];
set(axj, 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none')
set(axn, 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none')
hold on

% current flux tubes
n = 1;
in = false(size(j1));
out = false(size(j1));
flux_string_in = cell(1, ntubes);
flux_string_out = cell(1, ntubes);
k = 1;
for compare = double(do_compare):-1:0
    for tp = fieldnames(tube_pars)'
        tube_name = tp{1};
        pars = tube_pars.(tube_name);
        if compare
            tube = tubes_compare.(tube_name);
            color = 1 - pars.color;
        else
            tube = tubes.(tube_name);
            color = pars.color;
        end
        verts = tube.vertices;
        verts_2d = cellfun(@(v) v.*[0, 1, 1] + [x(end) - offs*k, 0, 0], verts, ...
            'UniformOutput', 0);
        points = tube.points;
        c0 = tube.caps.start;
        c1 = tube.caps.end;
        if tube.flux.in_axis == 3
            in = in | tube.flux.area.in;
        end
        if tube.flux.out_axis == 3
            out = out | tube.flux.area.out;
        end
        flux0 = tube.flux.in*scl.f;
        flux1 = tube.flux.out*scl.f;
        outline = tube.outline;
        inline = tube.inline;
        fprintf('Tube %s: influx = %.2f %s, outflux = %.2f %s, ratio = %.2f.\n' ...
            , tube_name, flux0, unt.f, flux1, unt.f, flux0 / flux1)
        flux_string_in{k} = ['\color[rgb]{', num2str(color), '}'...
            , num2str(flux0, '%3.2f'), '  ', '\color[rgb]{', num2str(cltx),'}'];
        flux_string_out{k} = ['\color[rgb]{', num2str(color), '}'...
            , num2str(flux1, '%3.2f'), '  ', '\color[rgb]{', num2str(cltx),'}'];
        k = k + 1;
        
        if pars.do_reverse
            clc0_tmp = clc1;
            clc1_tmp = clc0;
        else
            clc0_tmp = clc0;
            clc1_tmp = clc1;
        end
        if compare
            clc0_tmp = 1 - clc0_tmp;
            clc1_tmp = 1 - clc1_tmp;
        end

        stl = streamline(verts);
        if vind == 1 && pars.do_projection
            stl_2d = streamline(verts_2d);
        end
        plot3(c0(:, 1), c0(:, 2), c0(:, 3), 'Color', clc0_tmp);
        if vind == 1
            plot3(c0(:, 1), c0(:, 2), c0(:, 3)*0+z(1)+offs, 'Color', clc0_tmp, 'LineStyle', ':');
        end
        if ~iscell(c1)
            plot3(c1(:, 1), c1(:, 2), c1(:, 3), 'Color', clc1_tmp);
            if vind == 1
                plot3(c1(:, 1), c1(:, 2), c1(:, 3)*0+z(1)+offs, 'Color', clc1_tmp, 'LineStyle', ':');
            end
        else
            for i=1:length(c1)
                c = cell2mat(c1(i));
                plot3(c(:, 1), c(:, 2), c(:, 3), 'Color', clc1_tmp);
                if vind == 1
                    plot3(c(:, 1), c(:, 2), c(:, 3)*0+z(1), 'Color', clc1_tmp, 'LineStyle', ':');
                end
            end
        end
        set(stl, 'Color', [color, stlo], 'LineWidth', linw/2)
        if vind == 1 && pars.do_projection
            set(stl_2d, 'Color', [color, stlo/4], 'LineWidth', linw/2)
        end
    end
end

% field-aligned current slice
if vind ~=2
    j1_tmp = j1;
    j1_tmp(in | out) = j1(in | out) / 2;
    j1_slice = repmat(j1_tmp, [1, 1, length(z)]);
    slice(axj, Xm, Ym, Zm, -j1_slice, [], [], z(1))
    colormap(axj, colorcet(clm.j))
    clim(axj, lim.j)
    if vind ==1
        clb = colorbar(axj);
        clb.Color = cltx;
        clb.Label.String = sprintf('j_{||} (%s)', unt.j);
        clb.Label.Color = cltx;
        clb.Position = [clbx, clbh+2*(1-2*clbh)/3, 0.015, clbh];
    end
end

% electron density slice
if vind ~= 3
    ne_slice = permute(repmat(ne, [1, 1, length(x)]), [1, 3, 2]);
    slice(axn, Xm, Ym, Zm, ne_slice, x(end), [], [])
    colormap(axn, colorcet(clm.n))
    clim(axn, lim.n)
    if vind == 1
        clb = colorbar(axn);
        clb.Color = cltx;
        clb.Label.String = sprintf('log_{10} n_e (%s)', unt.n);
        clb.Label.Color = cltx;
        clb.Position = [clbx, (1-2*clbh)/3, 0.015, clbh];
    end
end

% track data
if vind ~= 2
    for s = sats
        x_tmp = track_x.(s);
        y_tmp = track_y.(s);
        z_tmp = ones(size(x_tmp)) * z(1);
        if strcmp(track_data_type, 'flow')
            vx_tmp = track_vx.(s) * scl.qv;
            vy_tmp = track_vy.(s) * scl.qv;
            vz_tmp = zeros(size(x_tmp));
        elseif strcmp(track_data_type, 'current')
            vx_tmp = track_fac.(s) * scl.qv;
            vy_tmp = zeros(size(x_tmp));
            vz_tmp = zeros(size(x_tmp));
        end
        quiver3(x_tmp, y_tmp, z_tmp, vx_tmp, vy_tmp, vz_tmp, '.-', 'Color', cltd)
    end
end

% electric field
if vind ~= 2
    n_qx = 3;
    n_qy = 4;
    qx = linspace(-(n_qx - 1) / 2, (n_qx - 1) / 2, n_qx) * range(x) / n_qx;
    qy = linspace(-(n_qy - 1) / 2, (n_qy - 1) / 2, n_qy) * range(y) / n_qy;
    [~, qx_ids] = min(abs(x - qx));
    [~, qy_ids] = min(abs(y - qy));
    Xm_q = Xm(qy_ids, qx_ids, 1:2);
    Ym_q = Ym(qy_ids, qx_ids, 1:2);
    Zm_q = Zm(qy_ids, qx_ids, 1:2);
    E2_q = E2(qy_ids, qx_ids, :) * scl.qe;
    E3_q = E3(qy_ids, qx_ids, :) * scl.qe;
    E1_q = zeros(size(E2_q)) * scl.qe;
    quiver3(Xm_q, Ym_q, Zm_q, E2_q, E3_q, E1_q, 0, 'Color', clef, 'MaxHeadSize', 1)
    quiver3(Xm_q(1, 1, 2), Ym_q(1, 1, 2), Zm_q(1, 1, 2), E2_bg * scl.qe, E3_bg * scl.qe, 0, 0, 'Color', cleb, 'MaxHeadSize', 1)
    if vind == idef
        % anne = vecnorm([E2_q(1, 1, 2), E3_q(1, 1, 2), E1_q(1, 1, 2)]);
        anne = vecnorm([E2_bg, E3_bg]);
        text(Xm_q(1, 1, 2)*1.1, Ym_q(1, 1, 2)*1.15, Zm_q(1, 1, 2)*1.02, ...
            sprintf('%.1f %s', anne, unt.e), 'Color', cleb, 'FontSize', fnts * 0.6, ...
            'HorizontalAlignment', 'right', 'BackgroundColor', [0, 0, 0, 0.5])
    end
end

% annotation
if vind == 1
    annotation('textbox', [0.03, 0.96, 0.01, 0.01], 'String', [...
        'In (', unt.f, '):  ', [flux_string_in{:}], newline, ...
        'Out (', unt.f, '):  ', [flux_string_out{:}], newline, ...
        ], 'FitBoxToText', 'on', 'EdgeColor', 'none', 'BackgroundColor', 'none', ...
        'FontSize', fnts, 'FontName', fntn, 'Color', cltx)
end

% adjust view positioning
view(axa, angl)
xlim(axa, 1.03*lim.x); ylim(axa, 1.03*lim.y); zlim(axa, lim.z.*[1, 1.01])
pbaspect(axj, ar); pbaspect(axn, ar)

xlabel(axj, sprintf('Mag. east (%s)', unt.x), 'Rotation', xrot)
ylabel(axj, sprintf('Mag. north (%s)', unt.x), 'Rotation', yrot)
zlabel(axj, sprintf('Mag. up (%s)', unt.x))
xlabel(axn, sprintf('Mag. east (%s)', unt.x), 'Rotation', xrot)
ylabel(axn, sprintf('Mag. north (%s)', unt.x), 'Rotation', yrot)
zlabel(axn, sprintf('Mag. up (%s)', unt.x))

camzoom(axj, zoom); camzoom(axn, zoom)
camdolly(axj, -panx, -pany, 0); camdolly(axn, -panx, -pany, 0)
box(axa, 'on')

set(fig, 'Color', clbg, 'InvertHardcopy', 'off')
set(axa, 'Color', 'none', 'GridColor', cltx, 'MinorGridColor', cltx, ...
    'XColor', cltx, 'YColor', cltx, 'ZColor', cltx)

% save figure
if save_plot
    [~, filename_tmp] = fileparts(direc);
    if ~isempty(sffx)
        filename = [filename_tmp, '_', char(sffx), '.png'];
    else
        filename = [filename_tmp, '.png'];
    end
    filename = fullfile(direc, 'plots3d', filename);
    fprintf('Saving %s\n', filename)
    print(fig, filename, '-dpng', '-r96');
    close all
end
end
end
