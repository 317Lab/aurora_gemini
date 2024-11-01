% Description:
%   Unpacks simulation directory in order to run aurogem.tools.replicate.
%
% Example usage:
%   aurogem.sim.replication('path-to-simulation')
%
% Arguments:
%   direc           simulation directory
%
% Contact:
%   jules.van.irsel.gr@dartmouth.edu
%
% Revisions:
%   07/23/2024  initial implementation (jvi)
%

function replication(direc)
arguments
    direc (1, :) char {mustBeFolder}
end

direc = fullfile(direc, 'ext');
cfg = gemini3d.read.config(direc);
xg = gemini3d.grid.cartesian(cfg);
x2 = xg.x2(3:end-2);

%% image
file_image = fullfile(direc, 'precipitation.h5');
try
    image.pos_type = h5readatt(file_image, '/', 'pos_type');
catch
    image.pos_type = 'angular';
end
if strcmp(image.pos_type, 'linear')
    pos1 = h5read(file_image, '/Coordinates/Magnetic/East');
    pos2 = h5read(file_image, '/Coordinates/Magnetic/North');
elseif strcmp(image.pos_type, 'angular')
    pos1 = h5read(file_image, '/Coordinates/Magnetic/Longitude');
    pos2 = h5read(file_image, '/Coordinates/Magnetic/Latitude');
end
[POS1, POS2] = ndgrid(pos1, pos2);

image.pos(:, :, 1) = POS1;
image.pos(:, :, 2) = POS2;
image.flux = h5read(file_image, '/Derived/Energy/Flux');
image.energy = h5read(file_image, '/Derived/Energy/Characteristic');
image.pedersen = h5read(file_image, '/Derived/Conductance/Pedersen');
image.hall = h5read(file_image, '/Derived/Conductance/Hall');

%% tracks
file_track = fullfile(direc, 'tracks.h5');
sats = char(cfg.used_tracks);
try
    tracks_pos_type = h5readatt(file_track, '/', 'pos_type');
catch
    tracks_pos_type = 'angular';
end

for sat = sats
    if strcmp(tracks_pos_type , 'linear')
        tmp.pos(:, 1) = h5read(file_track, ['/', sat, '/Coordinates/Magnetic/East']);
        tmp.pos(:, 2) = h5read(file_track, ['/', sat, '/Coordinates/Magnetic/North']);
    elseif strcmp(tracks_pos_type , 'angular')
        tmp.pos(:, 1) = h5read(file_track, ['/', sat, '/Coordinates/Magnetic/Longitude']);
        tmp.pos(:, 2) = h5read(file_track, ['/', sat, '/Coordinates/Magnetic/Latitude']);
    end
    tmp.flow(:, 1) = h5read(file_track, ['/', sat, '/Flow/Magnetic/East']);
    tmp.flow(:, 2) = h5read(file_track, ['/', sat, '/Flow/Magnetic/North']);
    tmp.fac(:, 1) = h5read(file_track, ['/', sat, '/Current/FieldAligned']);
    tmp.pos_type = tracks_pos_type;
    tracks.(sat) = tmp;
    clear('tmp')
end

%% replicate
plot_suffix = strrep(strip(char(cfg.plot_suffix), '_'), ' ', '');
if ~isempty(plot_suffix)
    plot_suffix = ['_', plot_suffix];
end
direc_out = fullfile(direc, 'output');
if ~exist(direc_out, 'dir')
    mkdir(direc_out)
end
if length(sats) == 2
    wsl = cfg.weighting_scale_length;
else
    wsl = 1;
end
try
    flow_bg = cfg.flow_background;
catch
    flow_bg = [nan, nan];
end
try
    swap_primary = cfg.swap_primary;
catch
    swap_primary = false;
end
try
    driver = cfg.driver;
catch
    driver = 'flow';
end

[bc, resnorm, E2_bg, E3_bg, v2_int, v3_int, weight0, bound] ...
    = aurogem.tools.replicate(tracks, image, xg ...
    , driver = driver ...
    , data_smoothing_window = cfg.data_smoothing_window ...
    , boundary_smoothing_window = cfg.boundary_smoothing_window ...
    , show_plots = true ...
    , save_plots = true ...
    , direc = direc_out ...
    , suffix = plot_suffix ...
    , add_phi_background = cfg.add_phi_background ...
    , fit_harmonic = cfg.fit_harmonic_function ...
    , num_replications = cfg.replication_number ...
    , arc_definition = cfg.arc_definition ...
    , edge_method = "contour" ...
    , do_rotate = true ...
    , do_scale = true ...
    , contour_values = cfg.contour_values ...
    , harmonic_mask = cfg.harmonic_mask ...
    , weighting_scale_length = wsl ...
    , flow_bg = flow_bg ...
    , swap_primary = swap_primary ...
    );

weight.(sats(1)) = weight0;
if length(sats) == 2
    weight.(sats(2)) = 1-weight0;
end
bound_prim = [x2; bound.A(x2)];
bound_scnd = [x2; bound.B(x2)];
h5make = @aurogem.tools.h5make;

answered = false;
while not(answered)
    answer = input('Save boundary condition? (Y/n) ', 's');
    if strcmp(answer, 'Y')
        answered = true;
        happy = true;
    elseif strcmp(answer, 'n')
        answered = true;
        happy = false;
    end
end

if happy
    if strcmp(driver, 'flow')
        file_bc = fullfile(direc, 'potential.h5');
    else
        file_bc = fullfile(direc, 'current.h5');
    end
    fprintf('Saving %s\n', file_bc)

    if strcmp(driver, 'flow')
        h5make(file_bc, '/Driver/Potential', bc, 'Convection electric potential' ...
            , units='Volts', size='lxp x lyp')
        h5make(file_bc, '/Driver/ResNorm', resnorm, 'Relative residual norm' ...
            , size='1 x 2 (East, North)', note='norm(residual)/norm(interpolated)')
        h5make(file_bc, '/InterpolatedFlow/East', v2_int, 'Magnetic eastward interpolated flow' ...
            , units='Meters/second', size='lxp x lyp')
        h5make(file_bc, '/InterpolatedFlow/North', v3_int, 'Magnetic northward interpolated flow' ...
            , units='Meters/second', size='lxp x lyp')
    else
        h5make(file_bc, '/Driver/Current', bc, 'Field aligned current' ...
            , units='Amperes/meter^2', size='lxp x lyp')
    end
    h5make(file_bc, '/BackgroundE/East', E2_bg ...
        , 'Magnetic eastward background electric field', units='Volts/meter')
    h5make(file_bc, '/BackgroundE/North', E3_bg ...
        , 'Magnetic northward background electric field', units='Volts/meter')
    h5make(file_bc, '/Boundary/Primary', bound_prim, 'Primary arc boundary'...
        , units='Meters', size='lxp x 2 (east, north)')
    h5make(file_bc, '/Boundary/Secondary', bound_scnd, 'Secondary arc boundary' ...
        , units='Meters', size='lxp x 2 (east, north)')
    for sat = sats
        h5make(file_bc, ['/Weights/', sat], weight.(sat), ['Track ', sat, ' weight map'] ...
            , size='lxp x lyp')
    end
end

close all
fclose all;
end