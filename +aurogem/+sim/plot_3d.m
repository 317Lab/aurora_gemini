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
%   04/26/2025  allow for start and end splits and kinks (jvi)
%

%% user input
sims = [...
"swop_20230210_35487_09_SD_AM_AC", ... #1
"swop_20230210_35487_09_PF_AM_AC", ... #2
"swop_20230210_35487_09_SD_UM_AC", ... #3
"swop_20230210_35487_09_SD_AM_xA", ... #4
...
"swop_20230212_37331_09_SD_AM_xC", ... #5
"swop_20230212_37331_09_PF_AM_xC", ... #6
"swop_20230212_37331_09_SD_UM_xC", ... #7
...
"swop_20230304_27012_09_SD_AM_xC", ... #8
"swop_20230304_27012_09_NB_AM_xC", ... #9
"swop_20230304_27012_09_SD_UM_xC", ... #10
...
"swop_20230304_36829_09_SD_AM_xB", ... #11
"swop_20230304_36829_09_SD_UM_xB", ... #12
...
"swop_20230314_24547_09_SD_AM_AC", ... #13
"swop_20230314_24547_09_PF_AM_AC", ... #14
"swop_20230314_24547_09_SD_AM_xA", ... #15
...
"swop_20230319_30210_09_SD_AM_xB", ... #16
"swop_20230319_30210_09_PF_AM_xB", ... #17
];

tube_sets = num2cell(ones(1, 17));
tube_sets{1} = [1, 2];
tube_sets{4} = 2;
num_tube_sets = numel([tube_sets{:}]);
tind = 0;

for sim_ind = 1
% clear('dat', 'xg')
sim = char(sims(sim_ind));
sim_compare = '';
% sim_compare = char(sims(4)); % when comparing tubes of two simulations

vinds = 1; % view angles to plot (1=iso, 2=side, 3=top)
spins = 0; %0:5:360; % for spinning animation (°)
save_plot = 0;
reload_tubes = true; % when only changing plotting parameters
reload_grid = false; % when only changing simulation data
reload_data = false; % when only changing tube configuration
publish_format = 'paper'; % publication format (paper, poster)
tube_set = tube_sets{sim_ind};
for tube_set_ind = tube_set
tic
tube_list = (1:3) + (tube_set_ind-1)*3; % list of flux tubes to plot
debug = 0;
do_contours = 0;
additional_sffx = num2str(tube_set_ind);

%% main
%#ok<*UNRCH>
for vind = vinds
for spin = spins
angl = [[-30, 32]; [-90, 0]; [0, 90]]; % view angle (°)
sffx = ["ISO", "SID", "TOP"];
fntn = 'Arial'; % font name
qntl = 0.95; % colorbar range quantile
clbh = 0.43; % colorbar height (relative)
stlo = 0.3 * (1 - debug); % flux tube opacity
offs = 0.5; % projection line offset (km)
idef = 3; % electric field legend vind (1, 3)
lnef = 20; % length of electric field background vector (km)
track_data_type = 'current'; % type of track data to plot ('current', 'flow')

clr.c0 = [0.0, 0.0, 0.0]; % start curve color (rgb)
clr.c1 = [0.0, 0.0, 0.0]; % end curve color (rgb)
clr.ef = [0.0, 0.0, 0.0]; % electric field color (rgb)
clr.eb = [1.0, 1.0, 0.0]; % electric field background (rgb)
clr.td = [1.0, 0.4, 1.0]; % track data color (rgb)

if strcmp(publish_format, 'paper')
    resf = 4;
    fnts = 10 * resf; % font size
    linw = 2; % line width
    pprw = [6.5, 2.58, 3.92] * resf; % paper width (inches)
    pprh = [5, 3.02, 3.02] * resf; % paper height (inches)
    clbx = 0.9; % colorbar horizontal position
    xrot = [18, 0, 0]; % x label rotation (deg)
    yrot = [-45, 0, 90]; % y label rotation (deg)
    sffx = sffx + "_P";

    clr.bg = [255, 255, 255] / 255; % background color (rgb)
    % clr.bg = [221, 221, 221] / 255;
    % clr.bg = [54, 54, 54] / 255;
    clr.tx = [0, 0, 0]; % text color (rgb)
    % clr.tx = [1, 1, 1];
else
    fnts = 18 * resf; % font size
    linw = 2; % line width
    pprw = [8.5, 4, 4] * resf; % paper width (inches)
    pprh = [7, 4, 3] * resf; % paper height (inches)
    clbx = 0.87; % colorbar horizontal position
    xrot = [18, 0, 0]; % x label rotation (deg)
    yrot = [-45, 0, 90]; % y label rotation (deg)

    clr.bg = [20, 21, 20] / 255; % background color (rgb)
    clr.tx = [1, 1, 1]; % text color (rgb)
end

% simulations directories
direc_root = getenv('GEMINI_SIM_ROOT');
direc = fullfile(direc_root, sim);
direc_compare = fullfile(direc_root, sim_compare);

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
    sffx = [sffx_tmp, '_COMPARE'];
elseif exist('spin', 'var')
    sffx = [sffx_tmp, '_', sprintf('%03i', spin)];
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
sats = strrep(sats{8}, 'x', '');

cfg = gemini3d.read.config(direc);
if not(exist('xg', 'var')) || reload_grid
    % xg = gemini3d.read.grid(direc);
    % xg = aurogem.tools.shrink(xg);
    xg = aurogem.grid.read(direc);
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
Ebg_norm = vecnorm([E2_bg, E3_bg]);
if Ebg_norm ~= 0
    scl.qe = lnef / Ebg_norm;
else
    scl.qe = 4;
end

% generate current flux tubes
if reload_tubes || not(exist('tubes', 'var'))
    for tp = fieldnames(tube_pars)'
        tube_name = tp{1};
        fprintf('Loading tube %s...\n', tube_name)
        tubes.(tube_name) = aurogem.tools.current_flux_tube(xg, dat, ...
            tube_pars.(tube_name), xlims = lim.x, ylims = lim.y, zlims = lim.z, debug=debug);
        if do_compare
            tubes_compare.(tube_name) = aurogem.tools.current_flux_tube(xg, dat_compare, ...
                tube_pars.(tube_name), xlims = lim.x, ylims = lim.y, zlims = lim.z, debug=debug);
        end
    end
end

if debug
    [~, xid] = min(abs(x-median(tubes.(tp{1}).caps.start(:, 1))));
    j2 = squeeze(dat.J2(lbz:ubz, xid, lby:uby))'*scl.j; % A/km^2 = uA/m^2
    [~, peak_z_id] = max(max(abs(j2)));
    [~, peak_y_id] = max(abs(j2(:, peak_z_id)));
    peak_j2 = j2(peak_y_id, peak_z_id);
    fprintf('Peak j2 of %.0f %s at (x, y, z) = (%.0f, %.0f, %.0f) %s\n', ...
        peak_j2, unt.j, x(xid), y(peak_y_id), z(peak_z_id), unt.x)
    ne = j2;
    lim.n = [-1, 1]*abs(peak_j2);
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
        fprintf('Tube %s: influx = %.2f %s, outflux = %.2f %s, ratio = %.2f.\n' ...
            , tube_name, flux0, unt.f, flux1, unt.f, flux0 / flux1)
        flux_string_in{k} = ['\color[rgb]{', num2str(color), '}'...
            , num2str(flux0, '%3.2f'), '  ', '\color[rgb]{', num2str(clr.tx),'}'];
        flux_string_out{k} = ['\color[rgb]{', num2str(color), '}'...
            , num2str(flux1, '%3.2f'), '  ', '\color[rgb]{', num2str(clr.tx),'}'];
        k = k + 1;

        stl = streamline(verts);
        if vind == 1 && pars.do_projection
            stl_2d = streamline(verts_2d);
        end
        cs = {c0, c1};
        for i = 1:2
            c = cs{i};
            if i == 1
                clr.tmp = clr.c0;
                lnwc = linw * 2;
            else
                clr.tmp = clr.c1;
                lnwc = linw;
            end
            if compare
                clr.tmp = 1-clr.tmp;
            end
            if ~iscell(c)
                plot3(c(:, 1), c(:, 2), c(:, 3), 'Color', clr.tmp, 'LineWidth', lnwc);
                if vind == 1
                    plot3(c(:, 1), c(:, 2), c(:, 3)*0+z(1)+offs, 'Color', clr.tmp, 'LineStyle', ':');
                end
            else
                for j=1:length(c)
                    cc = cell2mat(c(j));
                    plot3(cc(:, 1), cc(:, 2), cc(:, 3), 'Color', clr.tmp, 'LineWidth', lnwc);
                    if vind == 1
                        plot3(cc(:, 1), cc(:, 2), cc(:, 3)*0+z(1), 'Color', clr.tmp, 'LineStyle', ':');
                    end
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
if vind ~= 2
    j1_tmp = j1;
    j1_tmp(in | out) = j1(in | out) / 2;
    j1_slice = repmat(j1_tmp, [1, 1, length(z)]);
    slice(axj, Xm, Ym, Zm, -j1_slice, [], [], z(1))
    colormap(axj, colorcet(clm.j))
    clim(axj, lim.j)
    if vind ==1
        clb = colorbar(axj);
        clb.Color = clr.tx;
        clb.Label.String = sprintf('j_{||} (%s)', unt.j);
        clb.Label.Color = clr.tx;
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
        clb.Color = clr.tx;
        clb.Label.String = sprintf('log_{10} n_e (%s)', unt.n);
        clb.Label.Color = clr.tx;
        clb.Position = [clbx, (1-2*clbh)/3, 0.015, clbh];
    end
end

if do_contours > 0 && vind ~= 2
    % precip_file = fullfile(direc, cfg.prec_dir, ...
    %     gemini3d.datelab(cfg.times(end)) + ".h5");
    % Qp = h5read(precip_file, '/E0p');
    % Qp = Qp(lbx:ubx, lby:uby);
    [X, Y] = ndgrid(x, y);
    [~, zid] = min(abs(z - 110));
    ne_contour = log10(squeeze(dat.ne(zid, lbx:ubx, lby:uby))) * scl.n;
    iso_value = quantile(ne_contour(:), 0.7);
    ctrs = contour3(axn, X, Y, ne_contour, [1, 1] * 10.7);
    ctr_id = 1;
    while ctr_id < length(ctrs)
        ctr_len = ctrs(2, ctr_id);
        ctrs(:, ctr_id) = nan(2, 1);
        ctr_id = ctr_id + ctr_len + 1;
    end
    plot3(axn, ctrs(1, :), ctrs(2, :), ones(size(ctrs(1, :))) * min(z) * 1.1...
        , '--k')
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
        quiver3(x_tmp, y_tmp, z_tmp, vx_tmp, vy_tmp, vz_tmp, '.-', 'Color', clr.td)
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
    quiver3(Xm_q, Ym_q, Zm_q, E2_q, E3_q, E1_q, 0, 'Color', clr.ef, 'MaxHeadSize', 1)
    if Ebg_norm ~= 0
        quiver3(Xm_q(1, 1, 2), Ym_q(1, 1, 2), Zm_q(1, 1, 2) + 1, E2_bg * scl.qe, ...
            E3_bg * scl.qe, 0, 0, 'Color', clr.eb, 'MaxHeadSize', 1)
    end
    if vind == idef
        if Ebg_norm == 0
            e_anno = vecnorm([E2_q(1, 1, 2), E3_q(1, 1, 2)]) / scl.qe;
            cl.tmp = clr.ef;
        else
            e_anno = Ebg_norm;
            cl.tmp = clr.eb;
        end
        text(Xm_q(1, 1, 2)*1.05, Ym_q(1, 1, 2)*1.15, Zm_q(1, 1, 2)*1.02, ...
            sprintf('%.1f %s', e_anno, unt.e), 'Color', cl.tmp, 'FontSize', fnts * 0.6, ...
            'HorizontalAlignment', 'right', 'BackgroundColor', [[0, 0, 0] + ~Ebg_norm, 0.5])
    end
end

% annotation
if vind == 1
    annotation('textbox', [0.03, 0.96, 0.01, 0.01], 'String', [...
        'In (', unt.f, '):  ', [flux_string_in{:}], newline, ...
        'Out (', unt.f, '):  ', [flux_string_out{:}], newline, ...
        ], 'FitBoxToText', 'on', 'EdgeColor', 'none', 'BackgroundColor', 'none', ...
        'FontSize', fnts, 'FontName', fntn, 'Color', clr.tx)
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

set(fig, 'Color', clr.bg, 'InvertHardcopy', 'off')
set(axa, 'Color', 'none', 'GridColor', clr.tx, 'MinorGridColor', clr.tx, ...
    'XColor', clr.tx, 'YColor', clr.tx, 'ZColor', clr.tx)

% save figure
if save_plot
    [~, filename_tmp] = fileparts(direc);
    filename = [filename_tmp, '_', char(sffx), '.png'];
    if ~isempty(additional_sffx)
        filename = [filename(1:end-4), '_', additional_sffx, '.png'];
    end
    filename = fullfile(direc, 'plots3d', filename);
    fprintf('Saving %s\n', filename)
    print(fig, filename, '-dpng', '-r96');
    close all
end

end % spin
end % vind
tnow = toc;
if tind == 0
    tavg = tnow;
else
    tavg = (tavg * tind + tnow) / (tind + 1);
end
trem = (num_tube_sets - tind - 1) * tavg;
tind = tind + 1;
fprintf('------------------- sim_ind = %i, tube_set_ind = %i, %.1f minutes remaining -------------------\n', sim_ind, tube_set_ind, trem / 60)
end % tube_set_ind
end % sim_ind