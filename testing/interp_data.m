direc = fullfile('..', 'public_html', 'Gemini3D', 'swop_20230210_35487_A_07');
cfg = gemini3d.read.config(direc);
xg = gemini3d.read.grid(direc);
dat = gemini3d.read.frame(direc, 'time', cfg.times(end));

%%
x1 = double(xg.x1(3:end-2)) / 1e3;
x2 = double(xg.x2(3:end-2)) / 1e3;
x3 = double(xg.x3(3:end-2)) / 1e3;
[X1, X2, X3] = ndgrid(x1, x2, x3);
fv1 = griddedInterpolant(X1, X2, X3, dat.v1);
fv2 = griddedInterpolant(X1, X2, X3, dat.v2);
fv3 = griddedInterpolant(X1, X2, X3, dat.v3);

%%
direc_track = fullfile(direc, 'ext', 'tracks.h5');
track_x2 = h5read(direc_track, '/A/Coordinates/Magnetic/East') / 1e3;
track_x3 = h5read(direc_track, '/A/Coordinates/Magnetic/North') / 1e3;
track_x1 = 490 * ones(size(track_x2));

track_v1 = h5read(direc_track, '/A/Flow/Magnetic/Up')' / 1e3;
track_v2 = h5read(direc_track, '/A/Flow/Magnetic/East')' / 1e3;
track_v3 = h5read(direc_track, '/A/Flow/Magnetic/North')' / 1e3;

gem_v1 = fv1(track_x1, track_x2, track_x3) / 1e3;
gem_v2 = fv2(track_x1, track_x2, track_x3) / 1e3;
gem_v3 = fv3(track_x1, track_x2, track_x3) / 1e3;

%
close all
hold on
vs = 20;
vs1 = vs*10;
quiver3(track_x2, track_x3, track_x1, track_v2*vs, track_v3*vs, track_v1*vs, 0, '.-r')
quiver3(track_x2, track_x3, track_x1, gem_v2*vs1, gem_v3*vs1, gem_v1*vs1, 0, '.-b')
quiver3(min(x2)+50, min(x3)+30, track_x1(1), 1*vs, 0, 0, 0, '.-r')
quiver3(min(x2)+50, min(x3)+20, track_x1(1), 1*vs1, 0, 0, 0, '.-b')
text(min(x2)+50, min(x3)+40, '1 km/s')
xlim([min(x2), max(x2)])
ylim([min(x3), max(x3)])
zlim([min(x1), max(x1)])