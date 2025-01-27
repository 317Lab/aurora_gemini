% direc = fullfile('..', 'public_html', 'Gemini3D', 'swop_20230210_35487_A_09');
% direc = fullfile('..', 'public_html', 'Gemini3D', 'swop_20230212_37331_C_09');
direc = fullfile('..', 'public_html', 'Gemini3D', 'swop_20230314_24547_A_09');
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

[~, tmp] = fileparts(direc);
tmp = strsplit(tmp, '_');
sat = tmp{4};

%%
direc_track = fullfile(direc, 'ext', 'tracks.h5');
track_x2 = h5read(direc_track, sprintf('/%s/Coordinates/Magnetic/East', sat)) / 1e3;
track_x3 = h5read(direc_track, sprintf('/%s/Coordinates/Magnetic/North', sat)) / 1e3;
track_x1 = 490 * ones(size(track_x2));

track_v1 = h5read(direc_track, sprintf('/%s/Flow/Magnetic/Up', sat))' / 1e3;
track_v2 = h5read(direc_track, sprintf('/%s/Flow/Magnetic/East', sat))' / 1e3;
track_v3 = h5read(direc_track, sprintf('/%s/Flow/Magnetic/North', sat))' / 1e3;

gem_v1 = fv1(track_x1, track_x2, track_x3) / 1e3;
gem_v2 = fv2(track_x1, track_x2, track_x3) / 1e3;
gem_v3 = fv3(track_x1, track_x2, track_x3) / 1e3;

track_nan = track_x3 > max(x3) | track_x3 < min(x3);
track_x2(track_nan) = nan;
track_x3(track_nan) = nan;

lim.x1 = [min(x1), max(x1)];
lim.x2 = [min([x2; track_x2]), max([x2; track_x2])];
lim.x3 = [min(x3), max(x3)];

%
close all
hold on
vs = 30;
vs1 = vs;
quiver3(track_x2, track_x3, track_x1, track_v2*vs, track_v3*vs, track_v1*vs, 0, '.-r')
quiver3(track_x2, track_x3, track_x1, gem_v2*vs1, gem_v3*vs1, gem_v1*vs1, 0, '.-b')
quiver3(min(x2)+50, min(x3)+30, track_x1(1), 1*vs, 0, 0, 0, '.-r')
quiver3(min(x2)+50, min(x3)+20, track_x1(1), 1*vs1, 0, 0, 0, '.-b')
text(min(x2)+50, min(x3)+40, '1 km/s')
xlim(lim.x2)
ylim(lim.x3)
zlim(lim.x1)
xticks([min(x2), max(x2)])
yticks([min(x3), max(x3)])
grid on
legend('efi', 'gemini')
print(gcf, fullfile(direc, 'compare1.png'), '-dpng', '-r96')

%%
close all
hold on
plot(track_x3, track_v2, '-r')
plot(track_x3, track_v3, '-b')
plot(track_x3, gem_v2, '--r')
plot(track_x3, gem_v3, '--b')
legend('efi v_2', 'efi v_3', 'gemini v_2', 'gemini v_3')
ylim([-1,1] * max(abs([track_v1; track_v2; track_v3])) * 1.05)
xlim(lim.x3)
xticks([min(x3), max(x3)])
print(gcf, fullfile(direc, 'compare2.png'), '-dpng', '-r96')

close all