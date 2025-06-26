% sim = 'swop_20230210_35487_09_PF_AM_AC';
sim = 'swop_20230304_27012_09_NB_AM_xC';
direc_root = getenv('GEMINI_SIM_ROOT');
direc = fullfile(direc_root, sim);
direc_track = fullfile(direc, 'ext', 'tracks.h5');
direc_data = fullfile(direc, 'ext', 'data');

cfg = gemini3d.read.config(direc);
xg = aurogem.grid.read(direc);
dat = gemini3d.read.frame(direc, 'time', cfg.times(end), ...
    'vars', ["v2", "v3"]);

x1 = double(xg.x1(3:end-2))' / 1e3;
x2 = double(xg.x2(3:end-2))' / 1e3;
x3 = double(xg.x3(3:end-2))' / 1e3;
[X1, X2, X3] = ndgrid(x1, x2, x3);
fv2 = griddedInterpolant(X1, X2, X3, dat.v2);
fv3 = griddedInterpolant(X1, X2, X3, dat.v3);

%%
s = 'C';
track_x2 = h5read(direc_track, sprintf('/%s/Coordinates/Magnetic/East', s))' / 1e3;
track_x3 = h5read(direc_track, sprintf('/%s/Coordinates/Magnetic/North', s))' / 1e3;
track_x1 = 490 * ones(size(track_x2));

track_time = datetime(h5read(direc_track, sprintf('/%s/Time/Unix', s)), 'ConvertFrom', 'posixtime');
[~, track_id] = min(abs(track_time - cfg.times(end)));
n = (length(track_x1) - 1) / 2;
track_ids = track_id + (-n:n);
track_time = track_time(track_ids);

track_v2 = h5read(direc_track, sprintf('/%s/Flow/Magnetic/East', s));
track_v3 = h5read(direc_track, sprintf('/%s/Flow/Magnetic/North', s));

track_del_ids = track_x3 > max(x3) | track_x3 < min(x3);
track_time(track_del_ids) = [];
track_x1(track_del_ids) = [];
track_x2(track_del_ids) = [];
track_x3(track_del_ids) = [];
track_v2(track_del_ids) = [];
track_v3(track_del_ids) = [];

efi_cad = 2; % Hz
efi_vsat = 7628; % m/s
min_scale_length = 16e3; % m
efi_w = 2 * min_scale_length / (efi_vsat * efi_cad);
track_v2 = smoothdata(track_v2, 'gaussian', efi_w);
track_v3 = smoothdata(track_v3, 'gaussian', efi_w);

%%
direc_novx = fullfile(direc_data, dir([direc_data, filesep, sprintf('*EFI%s*novx.h5', s)]).name);
novx_time = datetime(h5read(direc_novx, '/Timestamp'), 'ConvertFrom', 'posixtime');
novx_v2 = h5read(direc_novx, '/ViMagE')';

novx_del_ids = novx_time < track_time(1) | novx_time > track_time(end);
novx_time(novx_del_ids) = [];
novx_v2(novx_del_ids) = [];
novx_v2 = smoothdata(novx_v2, 'gaussian', efi_w);
novx_x3 = interp1(track_time, track_x3, novx_time);

%%
gem_v2 = fv2(track_x1, track_x2, track_x3);
gem_v3 = fv3(track_x1, track_x2, track_x3);

close all
figure('PaperUnits', 'inches', 'PaperPosition', [0, 0, 6.5, 4] * 4)

hold on
plot(novx_x3, novx_v2, '-r')
plot(track_x3, gem_v2, '--r')
plot(track_x3, gem_v3, '--b')
legend('TII cross-track V_i', 'GEMINI v_2', 'GEMINI v_3')
ylim(gca, [-1, 1] * 5e2)
filename = '00_tii_example.png';
% print(gcf, fullfile('plots', filename), '-dpng', '-r96')

% close all