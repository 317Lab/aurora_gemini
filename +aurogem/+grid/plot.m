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
%   lon_lim = [0,360];      (option) longitude plotting limits
%   lat_lim = [41.7,90];    (option) latitude plotting limits
%   scale = 2               (option) plot scale
%   cam_orbit = [-60,-60]   (option) input to camorbit
%   cam_zoom = 3;           (option) input to camzoom
%   xg struct = struct      (option) simulation grid
%
% Contact:
%   jules.van.irsel.gr@dartmouth.edu
%
% Revisions:
%   07/23/2024  initial implementation (jvi)
%

function plot(direc,opts)
arguments
    direc (1,:) char {mustBeFolder}
    opts.plot_alaska (1,1) logical = true
    opts.plot_canada (1,1) logical = true
    opts.plot_cities (1,1) logical = true
    opts.plot_altitudes (1,1) logical = true
    opts.plot_xyz (1,1) logical = false
    opts.plot_earth (1,1) logical = false
    opts.lon_lim (1,2) double {mustBeNonnegative} = [0,360];
    opts.lat_lim (1,2) double {mustBeNonnegative} = [41.7,90];
    opts.scale (1,1) double {mustBePositive} = 2
    opts.cam_orbit (1,2) double = [-60,-60]
    opts.cam_zoom (1,1) double = 3;
    opts.xg struct = struct
end

%% parameters
scale = opts.scale;
RE = 6371;
offset = 80;
dart_green = [30,92,141]/255;

%% grid
if isempty(fieldnames(opts.xg))
    xg = gemini3d.read.grid(direc);
else
    xg = opts.xg;
end
az = deg2rad(xg.glon);
el = deg2rad(xg.glat);
ra = xg.alt/1e3 + RE;
[model_x,model_y,model_z] = sph2cart(az,el,ra);

%% earth
[earth_x,earth_y,earth_z] = sphere(50);
earth_x = RE*earth_x;
earth_y = RE*earth_y;
earth_z = RE*earth_z;

%% Alaska
lim.lon = opts.lon_lim;
lim.lat = opts.lat_lim;
if opts.plot_alaska
    shp.usa = shaperead('+aurogem/+grid/map_data/gadm41_USA_1.shp');
    [map_usa_x,map_usa_y,map_usa_z] = map(shp.usa,lim,RE);
end
if opts.plot_canada
    shp.cad = shaperead('+aurogem/+grid/map_data/gadm41_CAN_1.shp');
    [map_cad_x,map_cad_y,map_cad_z] = map(shp.cad,lim,RE);
end

%% points
mid_x = mean(model_x(:));
mid_y = mean(model_y(:));
mid_z = mean(model_z(:));

min_x = model_x(1,end,end);
min_y = model_y(1,end,end);
min_z = model_z(1,end,end);
max_x = model_x(end,end,end);
max_y = model_y(end,end,end);
max_z = model_z(end,end,end);

[fair_x,fair_y,fair_z] = sph2cart(deg2rad(-147.7200+360),deg2rad(64.8401),RE+10);
[anch_x,anch_y,anch_z] = sph2cart(deg2rad(-149.8997+360),deg2rad(61.2176),RE+10);

%% plot
close all
fts = 10*scale;
ftn = 'Arial';

reset(0)
set(0,'defaultFigurePaperUnits','inches')
aurogem.tools.setall(0,'FontName',ftn)
aurogem.tools.setall(0,'FontSize',fts)
aurogem.tools.setall(0,'Multiplier',1)
set(0,'defaultAxesFontSizeMode','manual')
set(0,'defaultSurfaceEdgeColor','flat')

figure
set(gcf,'PaperPosition',[0,0,6.5,6.5]*scale)
hold on
axis off
model(model_x,model_y,model_z,[0,0.7,0])
if opts.plot_earth
    surf(earth_x,earth_y,earth_z)
end
if opts.plot_alaska
    plot3(map_usa_x,map_usa_y,map_usa_z,'Color',[1,1,1]*0.4)
end
if opts.plot_canada
    plot3(map_cad_x,map_cad_y,map_cad_z,'Color',[1,1,1]*0.4)
end
if opts.plot_cities
    scatter3(fair_x,fair_y,fair_z,100,'filled','MarkerFaceColor',dart_green)
    scatter3(anch_x,anch_y,anch_z,100,'filled','MarkerFaceColor',dart_green)
    text(fair_x+offset,fair_y-offset,fair_z,'Fairbanks' ...
        ,'Color',dart_green,'FontSize',fts)
    text(anch_x+offset,anch_y-offset,anch_z,'Anchorage' ...
        ,'Color',dart_green,'FontSize',fts)
end
if opts.plot_altitudes
    scatter3(min_x,min_y,min_z,100,'filled','MarkerFaceColor','b')
    scatter3(max_x,max_y,max_z,100,'filled','MarkerFaceColor','b')
    text(min_x+offset,min_y-offset,min_z ...
        ,sprintf('%.0f km',min(xg.alt(:))/1e3),'Color','b','FontSize',fts)
    text(max_x+offset,max_y-offset,max_z ...
        ,sprintf('%.0f km',max(xg.alt(:))/1e3),'Color','b','FontSize',fts)
end
if opts.plot_xyz
    plot3([mid_x,mid_x+10*offset],[mid_y,mid_y],[mid_z,mid_z],'r')
    plot3([mid_x,mid_x],[mid_y,mid_y+10*offset],[mid_z,mid_z],'r')
    plot3([mid_x,mid_x],[mid_y,mid_y],[mid_z,mid_z+10*offset],'r')
    text(mid_x+11*offset,mid_y,mid_z,'X','Color','r','FontSize',fts)
    text(mid_x,mid_y+11*offset,mid_z,'Y','Color','r','FontSize',fts)
    text(mid_x,mid_y,mid_z+11*offset,'Z','Color','r','FontSize',fts)
end
pbaspect([1,1,1])
camtarget([mid_x,mid_y,mid_z])
camzoom(opts.cam_zoom)
camorbit(opts.cam_orbit(1),opts.cam_orbit(2))

    function model(x,y,z,c)
        lx1 = size(x,1);
        lx2 = size(x,2);
        lx3 = size(x,3);
        for i1 = [1,lx1]
            for i2 = [1,lx2]
                xx = squeeze(x(i1,i2,:));
                yy = squeeze(y(i1,i2,:));
                zz = squeeze(z(i1,i2,:));
                plot3(xx,yy,zz,'Color',c)
            end
        end
        for i1 = [1,lx1]
            for i3 = [1,lx3]
                xx = squeeze(x(i1,:,i3));
                yy = squeeze(y(i1,:,i3));
                zz = squeeze(z(i1,:,i3));
                plot3(xx,yy,zz,'Color',c)
            end
        end
        for i2 = [1,lx2]
            for i3 = [1,lx3]
                xx = squeeze(x(:,i2,i3));
                yy = squeeze(y(:,i2,i3));
                zz = squeeze(z(:,i2,i3));
                plot3(xx,yy,zz,'Color',c)
            end
        end
    end

    function [map_x,map_y,map_z] = map(shp,lim,radius)
        map_lon = [shp.X]'+360;
        map_lat = [shp.Y]';
        map_lon(map_lon > lim.lon(2)) = NaN;
        map_lon(map_lon < lim.lon(1)) = NaN;
        map_lat(map_lat > lim.lat(2)) = NaN;
        map_lat(map_lat < lim.lat(1)) = NaN;
        [map_x,map_y,map_z] = sph2cart(deg2rad(map_lon),deg2rad(map_lat),radius);
    end

end