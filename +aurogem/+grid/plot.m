% Description:
%   Generates gemini input fields and precipitation files from configuration and
%   grid structures. Used with aurogem.sim.replication outputs. Option
%   for "setup_functions" in "setup" namelist.
%
% Example usage:
%   aurogem.grid.plot(path-to-simulation)
%
% Arguments:
%   direc                   simulation directory
%   plot_alaska = true      (option) whether to plot Alaska
%   plot_canada = true      (option) whether to plot Canada
%   plot_cities = true      (option) whether to plot cities
%   plot_altitudes = true   (option) whether to plot min and max altitudes
%   plot_xyz = false        (option) whether to plot xyz triad
%   plot_earth = false      (option) whether to plot Earth
%   lon_lim = [0, 360];      (option) longitude plotting limits
%   lat_lim = [41.7, 90];    (option) latitude plotting limits
%   scale = 2               (option) plot scale
%   cam_orbit = [-60, -60]   (option) input to camorbit
%   cam_zoom = 3;           (option) input to camzoom
%   xg struct = struct      (option) simulation grid
%
% Contact:
%   jules.van.irsel.gr@dartmouth.edu
%
% Revisions:
%   07/23/2024  initial implementation (jvi)
%

function plot(direc, opts)
arguments
    direc (1, :) char {mustBeFolder}
    opts.plot_alaska (1, 1) logical = true
    opts.plot_canada (1, 1) logical = true
    opts.plot_cities (1, 1) logical = true
    opts.plot_altitudes (1, 1) logical = true
    opts.plot_xyz (1, 1) logical = false
    opts.plot_earth (1, 1) logical = false
    opts.lon_lim (1, 2) double {mustBeNonnegative} = [0, 360];
    opts.lat_lim (1, 2) double {mustBeNonnegative} = [41.7, 90];
    opts.scale (1, 1) double {mustBePositive} = 2
    opts.cam_orbit (1, 2) double = [-60, -60]
    opts.cam_zoom (1, 1) double = 8;
    opts.xg struct = struct
    opts.swarm_file (1, :) string
    opts.pfisr_vvels_file (1, :) char {mustBeFile}
    opts.save_plot (1, 1) logical = false
end

%% parameters
scale = opts.scale;
RE = 6371;
win = 250;
% color_cities = [1, 1, 1];
color_cities = [1, 1 ,1] * 0.4;
% color_background = [13, 53, 18] / 255;
color_background = [1.0, 1.0, 1.0];
% color_model = [0, 0.8, 1];
color_model = [0.0, 0.0, 0.0];
% color_top = [1.0, 0.7, 0.0];
color_top = [1.0, 0.0, 0.0];
% color_bottom = [71, 212, 90] / 255;
color_bottom = [0.0, 0.0, 1.0];
color_swarm = [1.0, 0.8, 0.0];
color_pfisr = [0.0, 0.8, 0.0];
color_sdarn = [1.0, 0.4, 0.0];

%% grid & config
cfg = gemini3d.read.config(direc);
if isempty(fieldnames(opts.xg))
    xg = gemini3d.read.grid(direc);
else
    xg = opts.xg;
end
az = deg2rad(xg.glon);
el = deg2rad(xg.glat);
ra = xg.alt / 1e3 + RE;
[model_x, model_y, model_z] = sph2cart(az, el, ra);

%% Earth
[earth_x, earth_y, earth_z] = sphere(50);
earth_x = RE*earth_x;
earth_y = RE*earth_y;
earth_z = RE*earth_z;

%% Alaska
lim.lon = opts.lon_lim;
lim.lat = opts.lat_lim;
if opts.plot_alaska
    shp.usa = shaperead('+aurogem/+grid/map_data/gadm41_USA_1.shp');
    [map_usa_x, map_usa_y, map_usa_z] = map(shp.usa, lim, RE);
end
if opts.plot_canada
    shp.cad = shaperead('+aurogem/+grid/map_data/gadm41_CAN_1.shp');
    [map_cad_x, map_cad_y, map_cad_z] = map(shp.cad, lim, RE);
end

%% points
mid_x = mean(model_x(:));
mid_y = mean(model_y(:));
mid_z = mean(model_z(:));

top_x = mean(model_x(end, :, :), 'all');
top_y = mean(model_y(end, :, :), 'all');
top_z = mean(model_z(end, :, :), 'all');
bottom_x = mean(model_x(1, :, :), 'all');
bottom_y = mean(model_y(1, :, :), 'all');
bottom_z = mean(model_z(1, :, :), 'all');

min_x = model_x(1, 1, end);
min_y = model_y(1, 1, end);
min_z = model_z(1, 1, end);
max_x = model_x(end, 1, end);
max_y = model_y(end, 1, end);
max_z = model_z(end, 1, end);

[fair_x, fair_y, fair_z] = sph2cart(deg2rad(-147.7200+360), deg2rad(64.8401), RE+10);
[anch_x, anch_y, anch_z] = sph2cart(deg2rad(-149.8997+360), deg2rad(61.2176), RE+10);

%% tracks
plot_swarm = false;
if isfield(opts, 'swarm_file')
    plot_swarm = true;
    t0 = mean(cfg.times);
    swarm_x = []; swarm_y = []; swarm_z = [];
    for sf = opts.swarm_file
        swarm_time = datetime(h5read(sf, '/Timestamp'), 'ConvertFrom', 'posixtime');
        [~, swarm_id] = min(abs(swarm_time - t0));
        swarm_ids = swarm_id + (-win:win);
        swarm_glon = h5read(sf, '/Longitude');
        swarm_glat = h5read(sf, '/GeodeticLatitude');
        swarm_galt = h5read(sf, '/GeodeticAltitude');
        swarm_az = deg2rad(swarm_glon(swarm_ids));
        swarm_el = deg2rad(swarm_glat(swarm_ids));
        swarm_ra = swarm_galt(swarm_ids) / 1e3 + RE;
        [swarm_x_tmp, swarm_y_tmp, swarm_z_tmp] = sph2cart(swarm_az, swarm_el, swarm_ra);
        swarm_x = [swarm_x, nan, swarm_x_tmp]; %#ok<*AGROW>
        swarm_y = [swarm_y, nan, swarm_y_tmp];
        swarm_z = [swarm_z, nan, swarm_z_tmp];
    end
end

plot_pfisr = false;
if isfield(opts, 'pfisr_vvels_file')
    plot_pfisr = true;
    pfisr_glon = h5read(opts.pfisr_vvels_file, '/VvelsGeoCoords/Longitude');
    pfisr_glat = h5read(opts.pfisr_vvels_file, '/VvelsGeoCoords/Latitude');
    pfisr_galt = h5read(opts.pfisr_vvels_file, '/VvelsGeoCoords/Altitude');
    [~, pfisr_id] = min(abs(pfisr_galt(1, :) - max(xg.alt(:))/1e3));
    pfisr_id = pfisr_id - 1;
    pfisr_az = deg2rad(pfisr_glon(:, pfisr_id));
    pfisr_el = deg2rad(pfisr_glat(:, pfisr_id));
    pfisr_ra = pfisr_galt(:, pfisr_id) + RE;
    [pfisr_x, pfisr_y, pfisr_z] = sph2cart(pfisr_az, pfisr_el, pfisr_ra);
end

%% plot
close all
fts = 7*scale;
ftn = 'Arial';
lnw = 1.5;

reset(0)
aurogem.tools.setall(0, 'FontName', ftn)
aurogem.tools.setall(0, 'FontSize', fts)
aurogem.tools.setall(0, 'Multiplier', 1)
set(0, 'defaultAxesFontSizeMode', 'manual')
set(0, 'defaultSurfaceEdgeColor', 'flat')

hold on
axis off
model_surf(model_x, model_y, model_z, color_model, color_bottom, color_top, lnw)
plot3([mid_x, mid_x], [mid_y, mid_y + 300], [mid_z + 50, mid_z + 100], 'Color', color_model, 'LineWidth', lnw, 'LineStyle', ':')
text(mid_x, mid_y + 670, mid_z + 120, 'Model space', 'Color', color_model)
plot3([top_x, top_x], [top_y - 50, top_y - 270], [top_z, top_z - 100], 'Color', color_top*0.9, 'LineWidth', lnw, 'LineStyle', ':')
text(top_x, top_y - 300, top_z - 100, 'Top-boundary', 'Color', color_top)
plot3([bottom_x, bottom_x], [bottom_y, bottom_y + 90], [bottom_z, bottom_z - 200], 'Color', color_bottom, 'LineWidth', lnw, 'LineStyle', ':')
text(bottom_x, bottom_y + 100, bottom_z - 230, 'DASC', 'Color', color_bottom, 'HorizontalAlignment', 'right')
if plot_swarm
    plot3(swarm_x, swarm_y, swarm_z, 'Color', color_swarm, 'LineWidth', lnw, 'LineStyle', '-')
    text(swarm_x(win), swarm_y(win) - 140, swarm_z(win) - 180, 'Swarm', 'Color', color_swarm)
end
if plot_pfisr
    plot3(pfisr_x, pfisr_y, pfisr_z, ...
        'Color', color_pfisr, 'LineWidth', lnw, 'LineStyle', '-')
    text(pfisr_x(end), pfisr_y(end) + 240, pfisr_z(end), 'PFISR', 'Color', color_pfisr)
end
for ii = [-1, 0, 1] * 50
    quiver3(-2800+abs(ii), -1900+abs(ii), 5600+250+ii, 0, 100, 50, 0, '.-', 'Color', color_sdarn, 'LineWidth', lnw, 'MarkerSize', 20, 'Marker', '.')
end
text(-2800+500, -1600, 5600, 'SuperDARN', 'Color', color_sdarn)
if opts.plot_earth
    surf(earth_x, earth_y, earth_z)
else
    surf(earth_x, earth_y, earth_z, 'FaceColor', 'none', 'EdgeColor', 'none')
end
if opts.plot_alaska
    plot3(map_usa_x, map_usa_y, map_usa_z, 'Color', [1, 1, 1]*0.5)
end
if opts.plot_canada
    plot3(map_cad_x, map_cad_y, map_cad_z, 'Color', [1, 1, 1]*0.5)
end
if opts.plot_cities
    scatter3(fair_x, fair_y, fair_z + 10, 50, 'filled', 'MarkerFaceColor', color_cities)
    scatter3(anch_x, anch_y, anch_z + 10, 50, 'filled', 'MarkerFaceColor', color_cities)
    text(fair_x, fair_y + 160, fair_z + 20, 'Fairbanks' ...
        , 'Color', color_cities, 'FontSize', fts * 0.8 , 'HorizontalAlignment', 'center')
    text(anch_x, anch_y + 180, anch_z + 20, 'Anchorage' ...
        , 'Color', color_cities, 'FontSize', fts * 0.8 , 'HorizontalAlignment', 'center')
end
if opts.plot_altitudes
    plot3([min_x, min_x], [min_y + 10, min_y + 160], [min_z, min_z - 30], 'Color', color_model, 'LineWidth', lnw, 'LineStyle', ':')
    plot3([max_x, max_x], [max_y + 10, max_y + 160], [max_z, max_z - 30], 'Color', color_model, 'LineWidth', lnw, 'LineStyle', ':')
    text(min_x, min_y + 180, min_z - 40 ...
        , sprintf('%.0f km', min(xg.alt(:))/1e3), 'Color', color_model, 'FontSize', fts ...
        , 'HorizontalAlignment', 'right')
    text(max_x, max_y + 180, max_z - 40 ...
        , sprintf('%.0f km', max(xg.alt(:))/1e3), 'Color', color_model, 'FontSize', fts ...
        , 'HorizontalAlignment', 'right')
end
if opts.plot_xyz
    plot3([mid_x, mid_x + 1800], [mid_y, mid_y], [mid_z, mid_z], 'r')
    plot3([mid_x, mid_x], [mid_y, mid_y + 1800], [mid_z, mid_z], 'r')
    plot3([mid_x, mid_x], [mid_y, mid_y], [mid_z, mid_z + 1800], 'r')
    text(mid_x + 1980, mid_y, mid_z, 'X', 'Color', 'r', 'FontSize', fts)
    text(mid_x, mid_y + 1980, mid_z, 'Y', 'Color', 'r', 'FontSize', fts)
    text(mid_x, mid_y, mid_z + 1980, 'Z', 'Color', 'r', 'FontSize', fts)
end

pbaspect([1, 1, 1])
camtarget([mid_x, mid_y, mid_z])
camzoom(opts.cam_zoom)
camorbit(opts.cam_orbit(1), opts.cam_orbit(2))

set(gcf, 'Color', 'none');
set(gca, 'Color', 'none');

if opts.save_plot
    filename = fullfile(direc, 'grid_context.png');
    exportgraphics(gcf, filename, 'Resolution', 96 * 2 * scale, ...
        'BackgroundColor', color_background)
    close all

    im = imread(filename);
    crop = [150, 120, 0, 0] * scale;
    im_crop = im(1+crop(1):end-crop(2), 1+crop(3):end-crop(4), :);
    imwrite(im_crop, filename)
end

    function model_surf(x, y, z, c1, c2, c3, lw)
        for j = [1, size(x, 2)]
            xx = squeeze(x([1, end], j, [1, end]));
            yy = squeeze(y([1, end], j, [1, end]));
            zz = squeeze(z([1, end], j, [1, end]));
            surf(xx, yy, zz, 'EdgeColor', 1-c1, 'LineWidth', lw / 3, ...
                'FaceColor', c1, 'FaceAlpha', 0.5, 'EdgeAlpha', 1)
        end
        for k = [1, size(x, 3)]
            xx = squeeze(x([1, end], [1, end], k));
            yy = squeeze(y([1, end], [1, end], k));
            zz = squeeze(z([1, end], [1, end], k));
            surf(xx, yy, zz, 'EdgeColor', 1-c1, 'LineWidth', lw / 3, ...
                'FaceColor', c1, 'FaceAlpha', 0.5, 'EdgeAlpha', 1)
        end
        for i = [1, size(x, 1)]
            if i == 1
                cc = c2;
            else
                cc = c3;
            end
            xx = squeeze(x(i, [1, end], [1, end]));
            yy = squeeze(y(i, [1, end], [1, end]));
            zz = squeeze(z(i, [1, end], [1, end])) + 1;
            surf(xx, yy, zz, 'EdgeColor', cc, 'LineWidth', lw, 'FaceColor', cc, 'FaceAlpha', 0.5)
        end
    end

    function [map_x, map_y, map_z] = map(shp, lim, radius)
        map_lon = [shp.X]'+360;
        map_lat = [shp.Y]';
        map_lon(map_lon > lim.lon(2)) = NaN;
        map_lon(map_lon < lim.lon(1)) = NaN;
        map_lat(map_lat > lim.lat(2)) = NaN;
        map_lat(map_lat < lim.lat(1)) = NaN;
        [map_x, map_y, map_z] = sph2cart(deg2rad(map_lon), deg2rad(map_lat), radius);
    end

end