% Description:
%   Text
%
% Example usage:
%   text
%
% Arguments:
%   text
%
% Contact:
%   jules.van.irsel.gr@dartmouth.edu
%
% Revisions:
%   04/10/2025  initial implementation (jvi)
%
function setup(data_direc, event_id, swarm_ids, opts)
arguments
    data_direc (1, :) char {mustBeFile}
    event_id (1, 1) int32 {mustBeNonnegative}
    swarm_ids (1, :) char {mustBeNonempty}
    opts.sim_version (1, 1) int32 {mustBeNonnegative} = 9
    opts.suffix (1, :) char = 'testing';
    opts.min_scale_length (1, 1) double {mustBePositive} = 16e3; % m
    opts.track_window (1, 1) double {mustBePositive} = 60; % s
    opts.track_shift (1, 2) double = [0, 0]; % CAUTION shifts track range(x3) * track_shift north
    opts.max_flow (1, 1) double {mustBeNonnegative} = 3e3; % m/s
    opts.driver (1, 1) string {mustBeMember(opts.driver, ["current", "flow"])}= 'current';
    opts.boundary_sim (1, :) char = '';
    opts.do_acc (1, 1) logical = true;
    opts.check_consecutive_efi (1, 1) logical = false;
    opts.background_flow_source (1, 1) string {mustBeMember( ...
        opts.background_flow_source, ["superdarn", "pfisr", "none"] ...
        )} = 'superdarn'; % superdarn, pfisr, none
    opts.Ap_override (1, 1) double = -1; % negative = no override
end

%% generate config structs
sim.tdur = 60; % s
sim.dtout = 5; % s
sim.alt_scale = '5150, 4850, 250000, 30000'; % m
sim.dtprec = 5; % s
sim.dtE0 = 5; % s
sim.PhiWBG = 0.1; % mW/m^2
sim.W0BG = 1e3; % eV

rep.driver = char(opts.driver);
rep.data_smoothing_window = -1;
rep.plot_suffix = '_';
rep.add_phi_background = 0;
rep.fit_harmonic_function = 1;
rep.replication_number = 400;
rep.arc_definition = 'Hall';
rep.harmonic_mask = [30e3, 30e3, 40e3]; % m
rep.used_tracks = swarm_ids;
rep.weighting_scale_length = 50e3; % m
rep.track_shift = opts.track_shift;

if any(rep.track_shift ~= 0)
    input('One or more tracks are shifted.')
end

grd.base.alt_min = 80e3;
grd.base.alt_max = 500e3;
grd.base.alt_scale = 20e3;
grd.base.Bincl = 90;
grd.base.lxp = 300;
grd.base.lyp = 400;

%% init
assert(not(ispc), 'Please run in Linux environment.')
root_sim = getenv('GEMINI_SIM_ROOT');
assert(~isempty(root_sim), ...
    'Add environment variable GEMINI_SIM_ROOT directing to gemini simulations')
h5make = @aurogem.tools.h5make;

suffix = opts.suffix;
if not(isempty(suffix))
    if not(strcmp(suffix(1), '_'))
        suffix = ['_', suffix];
    end
end

if opts.Ap_override >= 0
    suffix = ['_Apor', suffix];
end

suffix = ['_', pad(swarm_ids, 2, 'left', 'x'), suffix];

if opts.do_acc
    suffix = ['_AM', suffix];
else
    suffix = ['_UM', suffix];
end

bfs = opts.background_flow_source;
if strcmp(bfs, 'superdarn')
    suffix = ['_SD', suffix];
elseif strcmp(bfs, 'pfisr')
    suffix = ['_PF', suffix];
elseif strcmp(bfs, 'none')
    suffix = ['_NB', suffix];
else
    error('background_flow_source not found.')
end

events = readlines(fullfile(data_direc, 'events.dat'));
events_ids = find(startsWith(events, "------------"));
event_titles = strsplit(events(events_ids(1) + 1));
event = struct;
for e = events(events_ids(1):events_ids(2))'
    tmp = strsplit(e);
    if str2double(tmp{1}) == event_id
        for t = 1:length(event_titles)
            event.(event_titles(t)) = tmp{t};
        end
        break
    end
end
assert(~isempty(fieldnames(event)), 'Event ID %i not found.', event_id)

if isempty(opts.boundary_sim)
    rep.boundary_directory = filesep;
else
    rep.boundary_directory = fullfile(root_sim, opts.boundary_sim);
end

%% main
time = datetime(event.DATETIME, 'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss''Z');
fprintf('\n%#30s\n\n', sprintf(pad(' %s %s ', 80, 'both', '#'), time, swarm_ids))

rep.contour_values = str2double(strsplit(event.CONTOURS, ','));
rep.boundary_smoothing_window = str2double(event.BNDSMTH);
if strcmp(bfs, 'pfisr')
    rep.flow_background = str2double(strsplit(event.PFISRBG, ','));
elseif strcmp(bfs, 'superdarn')
    rep.flow_background = str2double(strsplit(event.SDARNBG, ','));
elseif strcmp(bfs, 'none')
    rep.flow_background = [0, 0];
end

time.Format = 'uuuuMMdd';
sim.x2parms = strrep(event.X2PARMS, ',', ', ');
sim.x3parms = strrep(event.X3PARMS, ',', ', ');
sim.eq_dir = fullfile(root_sim, 'ics', sprintf('swop_%s_%05i', ...
    time, second(time, 'secondofday')));
sim.flagdirich = strcmp(rep.driver, 'flow');
sim.W0_char = str2double(event.SRCE);

cfg_xg = grd.base;
cfg_xg.times = time;
cfg_xg.glat = str2double(event.GLAT);
cfg_xg.glon = str2double(event.GLON);
cfg_xg.xdist = str2double(event.XDIST);
cfg_xg.ydist = str2double(event.YDIST);
xg = gemini3d.grid.cartesian(cfg_xg);
xg.GLAT = squeeze(xg.glat(end, :, :));
xg.GLON = squeeze(xg.glon(end, :, :));
[X2, X3] = ndgrid(xg.x2(3:end-2), xg.x3(3:end-2));
geod_to_x2 = scatteredInterpolant(xg.GLON(:), xg.GLAT(:), X2(:), 'natural');
geod_to_x3 = scatteredInterpolant(xg.GLON(:), xg.GLAT(:), X3(:), 'natural');

% create simulation directories
direc_sim = fullfile(root_sim, sprintf('swop_%s_%05i_%02i%s', ...
    time, second(time, 'secondofday'), opts.sim_version, suffix));
if exist(direc_sim, 'dir')
    fprintf('Simulation directory exists:\n  %s\nwith contents:\n  %s\n', ...
        direc_sim, ls(direc_sim))
    ok = input('Continue? (y/n)', 's');
    if not(strcmp(ok, 'y'))
        return
    end
end
direc_ext = fullfile(direc_sim, 'ext');
direc_data = fullfile(direc_sim, 'ext', 'data');
if not(exist(direc_data, 'dir'))
    mkdir(direc_data)
end

% copy dasc data
fprintf('Copying DASC data files ...\n')
direc_dasc = fullfile(data_direc, 'dasc');
pattern_dasc = fullfile(direc_dasc, sprintf('INV_PKR_%s_%05i*.h5', ...
    time, second(time, 'secondofday')));
fns_dasc = {dir(pattern_dasc).name};
for f = fns_dasc
    if exist(fullfile(direc_data, f{1}), 'file'); continue; end
    copyfile(fullfile(direc_dasc, f{1}), fullfile(direc_data, f{1}), 'f')
end

% copy pfisr data
fprintf('Copying PFISR data files ...\n')
direc_pfisr = fullfile(data_direc, 'pfisr');
pattern_vvels = fullfile(direc_pfisr, sprintf('%s*vvels_lat.h5', time));
fns_pfisr = {dir(pattern_vvels).name};
for f = fns_pfisr
    if exist(fullfile(direc_data, f{1}), 'file'); continue; end
    copyfile(fullfile(direc_pfisr, f{1}), fullfile(direc_data, f{1}), 'f')
end

%% prepare precipitation.h5
fprintf('Preparing image ...\n')
if opts.do_acc
    image_fn = sprintf('INV_PKR_%s_%05i_ACC.h5', time, second(time, 'secondofday'));
else
    image_fn = sprintf('INV_PKR_%s_%05i.h5', time, second(time, 'secondofday'));
end
image_fn = fullfile(direc_data, image_fn);
image_glat = h5read(image_fn, '/Coordinates/Latitude');
image_glon = h5read(image_fn, '/Coordinates/Longitude');
image_flux = h5read(image_fn, '/Derived/Energy/Flux');
image_energy = h5read(image_fn, '/Derived/Energy/Characteristic');
image_pedersen = h5read(image_fn, '/Derived/Conductance/Pedersen');
image_hall = h5read(image_fn, '/Derived/Conductance/Hall');
image_red = h5read(image_fn, '/Optical/Red');
image_green = h5read(image_fn, '/Optical/Green');
image_blue = h5read(image_fn, '/Optical/Blue');

% interpolate onto geodetic grid
fflux = scatteredInterpolant(image_glon(:), image_glat(:), image_flux(:), 'natural');
fenergy = scatteredInterpolant(image_glon(:), image_glat(:), image_energy(:), 'natural');
fpedersen = scatteredInterpolant(image_glon(:), image_glat(:), image_pedersen(:), 'natural');
fhall = scatteredInterpolant(image_glon(:), image_glat(:), image_hall(:), 'natural');
fred = scatteredInterpolant(image_glon(:), image_glat(:), image_red(:), 'natural');
fgreen = scatteredInterpolant(image_glon(:), image_glat(:), image_green(:), 'natural');
fblue = scatteredInterpolant(image_glon(:), image_glat(:), image_blue(:), 'natural');

% smooth for minimum scale length
wx2 = 2 * opts.min_scale_length / median(xg.dx2h);
wx3 = 1;
image.flux = smoothdata2(fflux(xg.GLON, xg.GLAT), 'gaussian', {wx2, wx3});
image.energy = smoothdata2(fenergy(xg.GLON, xg.GLAT), 'gaussian', {wx2, wx3});
image.pedersen = smoothdata2(fpedersen(xg.GLON, xg.GLAT), 'gaussian', {wx2, wx3});
image.hall = smoothdata2(fhall(xg.GLON, xg.GLAT), 'gaussian', {wx2, wx3});
image.red = smoothdata2(fred(xg.GLON, xg.GLAT), 'gaussian', {wx2, wx3});
image.green = smoothdata2(fgreen(xg.GLON, xg.GLAT), 'gaussian', {wx2, wx3});
image.blue = smoothdata2(fblue(xg.GLON, xg.GLAT), 'gaussian', {wx2, wx3});
image.east = xg.x2(3:end-2);
image.north = xg.x3(3:end-2);
image.pos(:, :, 1) = X2;
image.pos(:, :, 2) = X3;
image.pos_type = 'linear';
image.time = time;

image.energy(image.energy < 200) = 200; % for ionize_fang:fang2010_spectrum

% save data
warning('on', 'aurogem:h5found')
lastwarn('')
prec_fn = fullfile(direc_ext, 'precipitation.h5');
fprintf('Saving %s ...\n', prec_fn)
h5make(prec_fn, '/Time/Year', int16(year(image.time)), 'Year', type='int16')

% suppress repetitive warnings
[~, warn_id] = lastwarn;
if strcmp(warn_id, 'aurogem:h5found')
    warning('off', warn_id)
    fprintf('Supressing further warnings of ID "aurogem:h5found"\n')
end

h5make(prec_fn, '/Time/DOY', int16(day(image.time, 'dayofyear')), 'Day of year', type='int16')
h5make(prec_fn, '/Time/Seconds', second(image.time, 'secondofday'), 'Seconds since midnight')
h5make(prec_fn, '/Time/Unix', posixtime(image.time), 'Unix time')

h5make(prec_fn, '/Coordinates/Magnetic/East', image.east, 'Magnetic eastward distance' ...
    , units='Meters', size='1 x Nx')
h5make(prec_fn, '/Coordinates/Magnetic/North', image.north, 'Magnetic northward distance' ...
    , units='Meters', size='1 x Ny')
h5make(prec_fn, '/Coordinates/Geodetic/Longitude', xg.GLON, 'Geodetic longitude' ...
    , units='Degrees east (0, 360)', size='Nx x Ny', foot_alt='110 km')
h5make(prec_fn, '/Coordinates/Geodetic/Latitude', xg.GLAT, 'Geodetic latitude' ...
    , units='Degrees north (-90, 90)', size='Nx x Ny', foot_alt='110 km')

h5make(prec_fn, '/Optical/Red', image.red, 'Red filter', units='Rayleighs', size='Nx x Ny')
h5make(prec_fn, '/Optical/Green', image.green, 'Green filter', units='Rayleighs', size='Nx x Ny')
h5make(prec_fn, '/Optical/Blue', image.blue, 'Blue filter', units='Rayleighs', size='Nx x Ny')

h5make(prec_fn, '/Derived/Energy/Flux', image.flux, 'Total precipitating energy flux' ...
    , units='Milliwatts/meter^2', size='Nx x Ny')
if opts.do_acc
    h5make(prec_fn, '/Derived/Energy/Characteristic', image.energy, 'Acceleration region potential drop' ...
        , units='Electronvolts', size='Nx x Ny')
else
    h5make(prec_fn, '/Derived/Energy/Characteristic', image.energy, 'Precipitating characteristic energy' ...
        , units='Electronvolts', size='Nx x Ny')
end
h5make(prec_fn, '/Derived/Conductance/Pedersen', image.pedersen, 'Pedersen conductance' ...
    , units='Siemens', size='Nx x Ny')
h5make(prec_fn, '/Derived/Conductance/Hall', image.hall, 'Hall conductance' ...
    , units='Siemens', size='Nx x Ny')

h5writeatt(prec_fn, '/', 'pos_type', 'linear')

%% prepare tracks.h5
fprintf('Copying Swarm data ...\n')
for sat = swarm_ids
    % copy swarm data
    direc_swarm = fullfile(data_direc, 'swarm');
    pattern_efi = fullfile(direc_swarm, sprintf('SW_*EFI%s*%s*.h5', sat, time));
    pattern_fac = fullfile(direc_swarm, sprintf('SW_*FAC%s*%s*.h5', sat, time));
    fns_swarm = {dir(pattern_efi).name, dir(pattern_fac).name};
    for f = fns_swarm
        if exist(fullfile(direc_data, f{1}), 'file'); continue; end
        copyfile(fullfile(direc_swarm, f{1}), fullfile(direc_data, f{1}), 'f')
    end

    for f = fns_swarm
        if contains(f{1}, 'FAC')
            fac_fn = fullfile(direc_data, f{1});
        elseif contains(f{1}, 'EXPT_EFI') && not(contains(f{1}, '_novx'))
            efi_fn = fullfile(direc_data, f{1});
        end
    end
    assert(contains([h5info(efi_fn).Datasets.Name], 'ViMagE'), 'Please run append_apex_data.py in swarm hdf5 data directory.')

    fprintf('Preparing tracks ...\n')

    fac_time = datetime(h5read(fac_fn, '/Timestamp'), 'ConvertFrom', 'posixtime');
    fac_cad = 1 / seconds(median(diff(fac_time)));
    [~, fac_id] = min(abs(fac_time - time));
    fac_ids = fac_id + (-opts.track_window * fac_cad : opts.track_window * fac_cad);
    fac_glat = h5read(fac_fn, '/GeodeticLatitude');
    fac_glon = wrapTo360(h5read(fac_fn, '/Longitude'));
    fac_glat110 = h5read(fac_fn, '/GeodeticLatitude110km');
    fac_glon110 = wrapTo360(h5read(fac_fn, '/Longitude110km'));
    fac_fac = h5read(fac_fn, '/FAC') * 1e-6;

    efi_time = datetime(h5read(efi_fn, '/Timestamp'), 'ConvertFrom', 'posixtime');
    efi_cad = 1 / seconds(median(diff(efi_time)));
    [~, efi_id] = min(abs(efi_time - time));
    efi_ids = efi_id + (-opts.track_window * efi_cad : opts.track_window * efi_cad);
    efi_gvu = h5read(efi_fn, '/ViU');
    efi_gve = h5read(efi_fn, '/ViE');
    efi_gvn = h5read(efi_fn, '/ViN');
    efi_gv = vecnorm([efi_gvu, efi_gve, efi_gvn]');
    efi_mvu = h5read(efi_fn, '/ViMagU');
    efi_mve = h5read(efi_fn, '/ViMagE');
    efi_mvn = h5read(efi_fn, '/ViMagN');
    efi_vsatu = -h5read(efi_fn, '/VsatC');
    efi_vsate = h5read(efi_fn, '/VsatE');
    efi_vsatn = h5read(efi_fn, '/VsatN');
    efi_vsat = vecnorm([efi_vsatu; efi_vsate; efi_vsatn]);

    % flow data processing
    efi_ids(efi_gv(efi_ids) > opts.max_flow) = [];
    if opts.check_consecutive_efi && (sum(diff(efi_ids) ~= 1) ~= 0)
        error('Non consecutive data')
    end

    % smooth for minimum scale length
    fac_w = 2 * opts.min_scale_length / median(efi_vsat(efi_ids) * fac_cad) / 2; % smooth half as much to balance against gradients
    efi_w = 2 * opts.min_scale_length / median(efi_vsat(efi_ids) * efi_cad);
    track.times = fac_time;
    track.glat = fac_glat(fac_ids);
    track.glon = fac_glon(fac_ids);
    track.glat110 = fac_glat110(fac_ids);
    track.glon110 = fac_glon110(fac_ids);
    track.east = smooth(geod_to_x2(track.glon110, track.glat110));
    track.north = smooth(geod_to_x3(track.glon110, track.glat110));
    track.fac = smoothdata(fac_fac(fac_ids), 'gaussian', fac_w);
    track.gvu = interp1(efi_time(efi_ids), ...
        smoothdata(efi_gvu(efi_ids), 'gaussian', efi_w), fac_time(fac_ids), 'spline');
    track.gve = interp1(efi_time(efi_ids), ...
        smoothdata(efi_gve(efi_ids), 'gaussian', efi_w), fac_time(fac_ids), 'spline');
    track.gvn = interp1(efi_time(efi_ids), ...
        smoothdata(efi_gvn(efi_ids), 'gaussian', efi_w), fac_time(fac_ids), 'spline');
    track.mvu = interp1(efi_time(efi_ids), ...
        smoothdata(efi_mvu(efi_ids), 'gaussian', efi_w), fac_time(fac_ids), 'spline');
    track.mve = interp1(efi_time(efi_ids), ...
        smoothdata(efi_mve(efi_ids), 'gaussian', efi_w), fac_time(fac_ids), 'spline');
    track.mvn = interp1(efi_time(efi_ids), ...
        smoothdata(efi_mvn(efi_ids), 'gaussian', efi_w), fac_time(fac_ids), 'spline');

    % save data
    swarm_fn = fullfile(direc_ext, 'tracks.h5');
    fprintf('Saving %s ...\n', swarm_fn)
    h5make(swarm_fn, ['/', sat, '/Time/Year'], int16(year(track.times)), 'Year', type='int16')
    h5make(swarm_fn, ['/', sat, '/Time/DOY'], int16(day(track.times, 'dayofyear')), 'Day of year', type='int16')
    h5make(swarm_fn, ['/', sat, '/Time/Seconds'], second(track.times, 'secondofday'), 'Seconds since midnight')
    h5make(swarm_fn, ['/', sat, '/Time/Unix'], posixtime(track.times), 'Unix time')

    h5make(swarm_fn, ['/', sat, '/Coordinates/Magnetic/East'], track.east, 'Magnetic eastward distance' ...
        , units='Meters')
    h5make(swarm_fn, ['/', sat, '/Coordinates/Magnetic/North'], track.north, 'Magnetic northward distance' ...
        , units='Meters')
    h5make(swarm_fn, ['/', sat, '/Coordinates/Geodetic/Longitude'], track.glon, 'Geodetic longitude' ...
        , units='Degrees east (0, 360)')
    h5make(swarm_fn, ['/', sat, '/Coordinates/Geodetic/Latitude'], track.glat, 'Geodetic latitude' ...
        , units='Degrees north (-90, 90)')
    h5make(swarm_fn, ['/', sat, '/Coordinates/Geodetic/FootLongitude'], track.glat110, 'Footpointed Geodetic longitude' ...
        , units='Degrees east (0, 360)', foot_alt='110 km')
    h5make(swarm_fn, ['/', sat, '/Coordinates/Geodetic/FootLatitude'], track.glat110, 'Footpointed Geodetic latitude' ...
        , units='Degrees north (-90, 90)', foot_alt='110 km')

    h5make(swarm_fn, ['/', sat, '/Flow/Magnetic/Up'], track.mvu, 'Magnetic upward plasma flow', units='Meters/second')
    h5make(swarm_fn, ['/', sat, '/Flow/Magnetic/East'], track.mve, 'Magnetic eastward plasma flow', units='Meters/second')
    h5make(swarm_fn, ['/', sat, '/Flow/Magnetic/North'], track.mvn, 'Magnetic northward plasma flow', units='Meters/second')
    h5make(swarm_fn, ['/', sat, '/Flow/Geodetic/Up'], track.gvu, 'Geodetic upward plasma flow', units='Meters/second')
    h5make(swarm_fn, ['/', sat, '/Flow/Geodetic/East'], track.gve, 'Geodetic eastward plasma flow', units='Meters/second')
    h5make(swarm_fn, ['/', sat, '/Flow/Geodetic/North'], track.gvn, 'Geodetic northward plasma flow', units='Meters/second')
    h5make(swarm_fn, ['/', sat, '/Current/FieldAligned'], track.fac, 'Field aligned current', units='Amperes/meter^2')
end

h5writeatt(swarm_fn, '/', 'pos_type', 'linear')

%% write main and replication config files
fprintf('Writing configuration files ...\n')
cfg_main = fullfile(direc_sim, 'config.nml');
cfg_rep = fullfile(direc_ext, 'config.nml');
fid0 = fopen(cfg_main, 'w');
fid1 = fopen(cfg_rep, 'w');
time.Format = 'uuuu,M,d';
[f107, ~, f107a, Ap] = aurogem.tools.activity(time);
if opts.Ap_override >= 0
    Ap = opts.Ap_override;
end

for fid = [fid0, fid1]
    fprintf(fid, '&base\n');
    fprintf(fid, sprintf('ymd = %s\n', time));
    fprintf(fid, sprintf('UTsec0 = %i\n', second(time, 'secondofday') - sim.tdur));
    fprintf(fid, sprintf('tdur = %i\n', sim.tdur));
    fprintf(fid, sprintf('dtout = %i\n', sim.dtout));
    fprintf(fid, sprintf('activ = %.1f,%.1f,%.0f\n', f107a, f107, Ap));
    fprintf(fid, 'tcfl = 0.9\n');
    fprintf(fid, 'Teinf = 1500\n');
    fprintf(fid, '/\n\n');

    fprintf(fid, '&setup\n');
    fprintf(fid, sprintf('glat = %.3f\n', cfg_xg.glat));
    fprintf(fid, sprintf('glon = %.2f\n', cfg_xg.glon));
    fprintf(fid, sprintf('xdist = %.0fe3\n', cfg_xg.xdist / 1e3));
    fprintf(fid, sprintf('ydist = %.0fe3\n', cfg_xg.ydist / 1e3));
    fprintf(fid, sprintf('alt_min = %.0fe3\n', cfg_xg.alt_min / 1e3));
    fprintf(fid, sprintf('alt_max = %.0fe3\n', cfg_xg.alt_max / 1e3));
    if fid==fid0
        fprintf(fid, sprintf('alt_scale = %s\n', sim.alt_scale));
        fprintf(fid, sprintf('x2parms = %s\n', sim.x2parms));
        fprintf(fid, sprintf('x3parms = %s\n', sim.x3parms));
        fprintf(fid, sprintf('lxp = %.0f\n', 1));
        fprintf(fid, sprintf('lyp = %.0f\n', 1));
        fprintf(fid, sprintf("eq_dir = '%s'\n", sim.eq_dir));
        fprintf(fid, "setup_functions = 'aurogem.functions.replicated'\n");
    else
        fprintf(fid, sprintf('alt_scale = %.0fe3\n', cfg_xg.alt_scale / 1e3));
        fprintf(fid, sprintf('lxp = %.0f\n', cfg_xg.lxp));
        fprintf(fid, sprintf('lyp = %.0f\n', cfg_xg.lyp));
    end
    fprintf(fid, sprintf('Bincl = %.0f\n', cfg_xg.Bincl));
    fprintf(fid, '/\n\n');

    fprintf(fid, '&flags\n');
    fprintf(fid, 'potsolve = 1\n');
    fprintf(fid, 'flagperiodic = 0\n');
    fprintf(fid, 'flagoutput = 1\n');
    fprintf(fid, '/\n\n');

    fprintf(fid, '&files\n');
    fprintf(fid, "file_format = 'h5'\n");
    fprintf(fid, "indat_size = 'inputs/simsize.h5'\n");
    fprintf(fid, "indat_grid = 'inputs/simgrid.h5'\n");
    fprintf(fid, "indat_file = 'inputs/initial_conditions.h5'\n");
    fprintf(fid, '/\n\n');

    if fid==fid0
        fprintf(fid, '&precip\n');
        fprintf(fid, sprintf('dtprec = %.0f\n', sim.dtprec));
        fprintf(fid, "prec_dir = 'inputs/particles'\n");
        fprintf(fid, '/\n\n');

        fprintf(fid, '&precip_BG\n');
        fprintf(fid, sprintf('PhiWBG = %.2f\n', sim.PhiWBG));
        fprintf(fid, sprintf('W0BG = %.0f\n', sim.W0BG));
        fprintf(fid, '/\n\n');

        fprintf(fid, '&efield\n');
        fprintf(fid, sprintf('dtE0 = %.0f\n', sim.dtE0));
        fprintf(fid, "E0_dir = 'inputs/fields'\n");
        fprintf(fid, '/\n\n');

        if opts.do_acc
            fprintf(fid, '&fang\n');
            fprintf(fid, 'flag_fang = 0\n');
            fprintf(fid, '/\n\n');

            fprintf(fid, '&fang_pars\n');
            fprintf(fid, 'diff_num_flux = 3\n');
            fprintf(fid, 'kappa = 1e4\n');
            fprintf(fid, 'bimax_frac = 1\n');
            fprintf(fid, sprintf('W0_char = %.0f\n', sim.W0_char));
            fprintf(fid, '/\n\n');
        end

        fprintf(fid, '&aurora_parameters\n');
        fprintf(fid, sprintf('flagdirich = %.0f\n', sim.flagdirich));
        fprintf(fid, '/\n');
    else
        fprintf(fid, '&replication\n');
        fprintf(fid, sprintf("driver = '%s'\n", rep.driver));
        fprintf(fid, sprintf("boundary_directory = '%s'\n", rep.boundary_directory));
        fprintf(fid, sprintf('data_smoothing_window = %i\n', rep.data_smoothing_window));
        fprintf(fid, sprintf('boundary_smoothing_window = %i\n', rep.boundary_smoothing_window));
        fprintf(fid, sprintf("plot_suffix = '%s'\n", rep.plot_suffix));
        fprintf(fid, sprintf('add_phi_background = %i\n', rep.add_phi_background));
        fprintf(fid, sprintf('fit_harmonic_function = %i\n', rep.fit_harmonic_function));
        fprintf(fid, sprintf('replication_number = %i\n', rep.replication_number));
        fprintf(fid, sprintf("arc_definition = '%s'\n", rep.arc_definition));
        fprintf(fid, sprintf('contour_values = %.1f,%.1f\n', rep.contour_values(1), rep.contour_values(2)));
        fprintf(fid, sprintf('harmonic_mask = %.0fe3,%.0fe3,%.0fe3\n', ...
            rep.harmonic_mask(1) / 1e3, rep.harmonic_mask(2) / 1e3, rep.harmonic_mask(3) / 1e3));
        fprintf(fid, sprintf("used_tracks = '%s'\n", rep.used_tracks));
        fprintf(fid, sprintf('flow_background = %.0f,%.0f\n', rep.flow_background(1), rep.flow_background(2)));
        fprintf(fid, sprintf('weighting_scale_length = %.0fe3\n', rep.weighting_scale_length / 1e3));
        fprintf(fid, sprintf('track_shift = %.2f,%.2f\n', rep.track_shift(1), rep.track_shift(2)));
        fprintf(fid, '/\n');
    end
end

fclose all;

%% run initial condition
status = aurogem.sim.run_ic(direc_sim);
if status ~= 0; return; end

%% perform replication
fprintf('Performing replication ...\n')
aurogem.sim.replication(direc_sim)

%% gemini setup
while true
    run_setup = input('Run GEMINI setup? (y/n) ', 's');
    if strcmp(run_setup, 'y')
        nas_username = input('NAS Agency User ID: ', 's');
        nas_group = input('NAS Group ID (GID): ', 's');
        slurm_account = input('sbatch account: ', 's');
        fprintf('Performing GEMINI setup ...\n')
        aurogem.sim.setup(direc_sim, nas_username, 'both', ...
            group_list=nas_group, slurm_account=slurm_account)
        break
    elseif strcmpi(run_setup, 'n')
        break
    end
end


end