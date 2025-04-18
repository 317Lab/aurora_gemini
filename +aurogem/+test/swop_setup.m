%%% OLD DONT USE
%% user input
direc = fullfile('..', 'Jules', 'thesis', 'data', 'paper2');
events_ids = 9;
sim_version = 9;
suffix = '_Ap';
boundary_sim = '';
do_acc = true;
do_bg_flow = true;
min_scale_length = 16e3; % m
track_window = 30; % s
max_flow = 5e3; % m/s
driver = 'current'; % flow or current

sim.tdur = 60; % s
sim.dtout = 5; % s
sim.alt_scale = '5150, 4850, 250000, 30000'; % m
sim.dtprec = 5; % s
sim.dtE0 = 5; % s
sim.PhiWBG = 0.1; % mW/m^2
sim.W0BG = 1e3; % eV

rep.driver = driver;
rep.data_smoothing_window = -1;
rep.plot_suffix = '_';
rep.add_phi_background = 0;
rep.fit_harmonic_function = 1;
rep.replication_number = 400;
rep.arc_definition = 'Hall';
rep.harmonic_mask = [30e3, 30e3, 40e3]; % m
rep.used_tracks = 'A';
%#ok<*UNRCH>

%% init
assert(ispc, 'Please run in Linux environment.')
if not(isempty(suffix))
    if not(strcmp(suffix(1), '_'))
        suffix = ['_', suffix];
    end
end
events = readlines(fullfile(direc, 'event_data.txt'));
events = events(2:end-1)';
root_sim = getenv('GEMINI_SIM_ROOT');
assert(~isempty(root_sim), ...
    'Add environment variable GEMINI_SIM_ROOT directing to gemini simulations')

h5make = @aurogem.tools.h5make;

grd.base.alt_min = 80e3;
grd.base.alt_max = 500e3;
grd.base.alt_scale = 20e3;
grd.base.Bincl = 90;
grd.base.lxp = 300;
grd.base.lyp = 400;

if isempty(boundary_sim)
    rep.boundary_directory = filesep;
else
    rep.boundary_directory = fullfile(root_sim, boundary_sim);
end

%% main
for e = events(events_ids)
    tmp = strsplit(e);
    time = datetime(tmp{2}, 'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss''Z');
    sat = char(tmp{3});
    fprintf('\n%#30s\n\n', sprintf(pad(' %s ', 80, 'both', '#'), time))
    
    rep.contour_values = str2double(strsplit(tmp{11}, ','));
    rep.flow_background = str2double(strsplit(tmp{12}, ',')) * do_bg_flow;
    rep.boundary_smoothing_window = str2double(tmp{14});
    
    time.Format = 'uuuuMMdd';
    sim.x2parms = strrep(tmp{8}, ',', ', ');
    sim.x3parms = strrep(tmp{9}, ',', ', ');
    sim.eq_dir = fullfile(root_sim, 'ics', sprintf('swop_%s_%05i', time, second(time, 'secondofday')));
    sim.flagdirich = strcmp(driver, 'flow');
    sim.W0_char = str2double(tmp{13});

    cfg_xg = grd.base;
    cfg_xg.times = time;
    cfg_xg.glat = str2double(tmp{4});
    cfg_xg.glon = str2double(tmp{5});
    cfg_xg.xdist = str2double(tmp{6});
    cfg_xg.ydist = str2double(tmp{7});
    xg = gemini3d.grid.cartesian(cfg_xg);
    xg.GLAT = squeeze(xg.glat(end, :, :));
    xg.GLON = squeeze(xg.glon(end, :, :));
    [X2, X3] = ndgrid(xg.x2(3:end-2), xg.x3(3:end-2));
    geod_to_x2 = scatteredInterpolant(xg.GLON(:), xg.GLAT(:), X2(:), 'natural');
    geod_to_x3 = scatteredInterpolant(xg.GLON(:), xg.GLAT(:), X3(:), 'natural');

    % create simulation directories
    direc_sim = fullfile(root_sim, sprintf('swop_%s_%05i_%s_%02i%s', ...
        time, second(time, 'secondofday'), sat, sim_version, suffix));
    direc_ext = fullfile(direc_sim, 'ext');
    direc_data = fullfile(direc_sim, 'ext', 'data');
    if not(exist(direc_data, 'dir'))
        mkdir(direc_data)
    end

    % copy dasc data
    direc_dasc = fullfile(direc, 'dasc_data');
    pattern_dasc = fullfile(direc_dasc, sprintf('INV_PKR_%s_%05i*.h5', time, second(time, 'secondofday')));
    fns_dasc = {dir(pattern_dasc).name};
    for f = fns_dasc
        if exist(fullfile(direc_data, f{1}), 'file'); continue; end
        copyfile(fullfile(direc_dasc, f{1}), fullfile(direc_data, f{1}), 'f')
    end

    % copy swarm data
    direc_swarm = fullfile(direc, 'swarm_data');
    pattern_efi = fullfile(direc_swarm, sprintf('SW_*EFI%s*%s*.h5', sat, time));
    pattern_fac = fullfile(direc_swarm, sprintf('SW_*FAC%s*%s*.h5', sat, time));
    fns_swarm = {dir(pattern_efi).name, dir(pattern_fac).name};
    for f = fns_swarm
        if exist(fullfile(direc_data, f{1}), 'file'); continue; end
        copyfile(fullfile(direc_swarm, f{1}), fullfile(direc_data, f{1}), 'f')
    end

    % copy pfisr data
    direc_pfisr = fullfile(direc, 'pfisr_data');
    pattern_vvels = fullfile(direc_pfisr, sprintf('%s*vvels_lat.h5', time));
    fns_pfisr = {dir(pattern_vvels).name};
    for f = fns_pfisr
        if exist(fullfile(direc_data, f{1}), 'file'); continue; end
        copyfile(fullfile(direc_pfisr, f{1}), fullfile(direc_data, f{1}), 'f')
    end

    %% prepare precipitation.h5
    if do_acc
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
    wx2 = 2 * min_scale_length / median(xg.dx2h);
    % wx3 = 2 * min_scale_length / median(xg.dx3h);
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
    prec_fn = fullfile(direc_ext, 'precipitation.h5');
    h5make(prec_fn, '/Time/Year', year(image.time), 'Year', type='int16')
    h5make(prec_fn, '/Time/DOY', day(image.time, 'dayofyear'), 'Day of year', type='int16')
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
    if do_acc
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
    for f = fns_swarm
        if contains(f{1}, 'FAC')
            fac_fn = fullfile(direc_data, f{1});
        elseif contains(f{1}, 'EXPT_EFI') && not(contains(f{1}, '_novx'))
            efi_fn = fullfile(direc_data, f{1});
        end
    end

    fac_time = datetime(h5read(fac_fn, '/Timestamp'), 'ConvertFrom', 'posixtime');
    fac_cad = 1 / seconds(median(diff(fac_time)));
    [~, fac_id] = min(abs(fac_time - time));
    fac_ids = fac_id + (-track_window * fac_cad : track_window * fac_cad);
    fac_glat = h5read(fac_fn, '/GeodeticLatitude');
    fac_glon = wrapTo360(h5read(fac_fn, '/Longitude'));
    fac_glat110 = h5read(fac_fn, '/GeodeticLatitude110km');
    fac_glon110 = wrapTo360(h5read(fac_fn, '/Longitude110km'));
    fac_mlat = h5read(fac_fn, '/MagneticLatitude');
    fac_mlon = wrapTo360(h5read(fac_fn, '/MagneticLongitude'));
    fac_fac = h5read(fac_fn, '/FAC') * 1e-6;

    efi_time = datetime(h5read(efi_fn, '/Timestamp'), 'ConvertFrom', 'posixtime');
    efi_cad = 1 / seconds(median(diff(efi_time)));
    [~, efi_id] = min(abs(efi_time - time));
    efi_ids = efi_id + (-track_window * efi_cad : track_window * efi_cad);
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
    efi_ids(efi_gv(efi_ids) > max_flow) = [];
    if sum(diff(efi_ids) ~= 1) ~= 0
        error('Non consecutive data')
    end
    
    % smooth for minimum scale length
    fac_w = 2 * min_scale_length / median(efi_vsat(efi_ids) * fac_cad) / 2; % smooth half as much to balance against gradients
    efi_w = 2 * min_scale_length / median(efi_vsat(efi_ids) * efi_cad);
    track.times = fac_time;
    track.glat = fac_glat(fac_ids);
    track.glon = fac_glon(fac_ids);
    track.glat110 = fac_glat110(fac_ids);
    track.glon110 = fac_glon110(fac_ids);
    track.east = smooth(geod_to_x2(track.glon110, track.glat110));
    track.north = smooth(geod_to_x3(track.glon110, track.glat110));
    track.fac = smoothdata(fac_fac(fac_ids), 'gaussian', fac_w);
    track.gvu = interp1(efi_time(efi_ids), ...
        smoothdata(efi_gvu(efi_ids), 'gaussian', efi_w), fac_time(fac_ids), 'cubic');
    track.gve = interp1(efi_time(efi_ids), ...
        smoothdata(efi_gve(efi_ids), 'gaussian', efi_w), fac_time(fac_ids), 'cubic');
    track.gvn = interp1(efi_time(efi_ids), ...
        smoothdata(efi_gvn(efi_ids), 'gaussian', efi_w), fac_time(fac_ids), 'cubic');
    track.mvu = interp1(efi_time(efi_ids), ...
        smoothdata(efi_mvu(efi_ids), 'gaussian', efi_w), fac_time(fac_ids), 'cubic');
    track.mve = interp1(efi_time(efi_ids), ...
        smoothdata(efi_mve(efi_ids), 'gaussian', efi_w), fac_time(fac_ids), 'cubic');
    track.mvn = interp1(efi_time(efi_ids), ...
        smoothdata(efi_mvn(efi_ids), 'gaussian', efi_w), fac_time(fac_ids), 'cubic');

    % save data
    flow_fn = fullfile(direc_ext, 'tracks.h5');
    h5make(flow_fn, ['/', sat, '/Time/Year'], year(track.times), 'Year', type='int16')
    h5make(flow_fn, ['/', sat, '/Time/DOY'], day(track.times, 'dayofyear'), 'Day of year', type='int16')
    h5make(flow_fn, ['/', sat, '/Time/Seconds'], second(track.times, 'secondofday'), 'Seconds since midnight')
    h5make(flow_fn, ['/', sat, '/Time/Unix'], posixtime(track.times), 'Unix time')
     
    h5make(flow_fn, ['/', sat, '/Coordinates/Magnetic/East'], track.east, 'Magnetic eastward distance' ...
        , units='Meters')
    h5make(flow_fn, ['/', sat, '/Coordinates/Magnetic/North'], track.north, 'Magnetic northward distance' ...
        , units='Meters')
    h5make(flow_fn, ['/', sat, '/Coordinates/Geodetic/Longitude'], track.glon, 'Geodetic longitude' ...
        , units='Degrees east (0, 360)')
    h5make(flow_fn, ['/', sat, '/Coordinates/Geodetic/Latitude'], track.glat, 'Geodetic latitude' ...
        , units='Degrees north (-90, 90)')
    h5make(flow_fn, ['/', sat, '/Coordinates/Geodetic/FootLongitude'], track.glat110, 'Footpointed Geodetic longitude' ...
        , units='Degrees east (0, 360)', foot_alt='110 km')
    h5make(flow_fn, ['/', sat, '/Coordinates/Geodetic/FootLatitude'], track.glat110, 'Footpointed Geodetic latitude' ...
        , units='Degrees north (-90, 90)', foot_alt='110 km')
    
    h5make(flow_fn, ['/', sat, '/Flow/Magnetic/Up'], track.mvu, 'Magnetic upward plasma flow', units='Meters/second')
    h5make(flow_fn, ['/', sat, '/Flow/Magnetic/East'], track.mve, 'Magnetic eastward plasma flow', units='Meters/second')
    h5make(flow_fn, ['/', sat, '/Flow/Magnetic/North'], track.mvn, 'Magnetic northward plasma flow', units='Meters/second')
    h5make(flow_fn, ['/', sat, '/Flow/Geodetic/Up'], track.gvu, 'Geodetic upward plasma flow', units='Meters/second')
    h5make(flow_fn, ['/', sat, '/Flow/Geodetic/East'], track.gve, 'Geodetic eastward plasma flow', units='Meters/second')
    h5make(flow_fn, ['/', sat, '/Flow/Geodetic/North'], track.gvn, 'Geodetic northward plasma flow', units='Meters/second')
    h5make(flow_fn, ['/', sat, '/Current/FieldAligned'], track.fac, 'Field aligned current', units='Amperes/meter^2')

    h5writeatt(flow_fn, '/', 'pos_type', 'linear')

    %% write main and replication config files
    cfg_main = fullfile(direc_sim, 'config.nml');
    cfg_rep = fullfile(direc_ext, 'config.nml');
    fid0 = fopen(cfg_main, 'w');
    fid1 = fopen(cfg_rep, 'w');
    time.Format = 'uuuu,M,d';
    [f107, ~, f107a, Ap] = aurogem.tools.activity(time);
    
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
            
            if do_acc
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
            fprintf(fid, sprintf('contour_values = %i,%i\n', rep.contour_values(1), rep.contour_values(2)));
            fprintf(fid, sprintf('harmonic_mask = %.0fe3,%.0fe3,%.0fe3\n', ...
                rep.harmonic_mask(1) / 1e3, rep.harmonic_mask(2) / 1e3, rep.harmonic_mask(3) / 1e3));
            fprintf(fid, sprintf("used_tracks = '%s'\n", rep.used_tracks));
            fprintf(fid, sprintf('flow_background = %.0f,%.0f\n', rep.flow_background(1), rep.flow_background(2)));
            fprintf(fid, '/\n');
        end
    end
end

fclose all;