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
    direc (1,:) char {mustBeFolder}
end

direc = fullfile(direc,'ext');
cfg = gemini3d.read.config(direc);
xg = gemini3d.grid.cartesian(cfg);
x2 = xg.x2(3:end-2);

%% image
file_image = fullfile(direc,'precipitation.h5');
mlon = h5read(file_image,'/Coordinates/Magnetic/Longitude');
mlat = h5read(file_image,'/Coordinates/Magnetic/Latitude');
[MLON,MLAT] = ndgrid(mlon,mlat);

image.pos(:,:,1) = MLON;
image.pos(:,:,2) = MLAT;
image.flux = h5read(file_image,'/Derived/Energy/Flux');
image.energy = h5read(file_image,'/Derived/Energy/Characteristic');
image.pedersen = h5read(file_image,'/Derived/Conductance/Pedersen');
image.hall = h5read(file_image,'/Derived/Conductance/Hall');

%% tracks
file_track = fullfile(direc,'tracks.h5');
sats = char(cfg.used_tracks);

for sat = sats
    tmp.pos(:,1) = h5read(file_track,['/',sat,'/Coordinates/Magnetic/Longitude']);
    tmp.pos(:,2) = h5read(file_track,['/',sat,'/Coordinates/Magnetic/Latitude']);
    tmp.flow(:,1) = h5read(file_track,['/',sat,'/Flow/Magnetic/East']);
    tmp.flow(:,2) = h5read(file_track,['/',sat,'/Flow/Magnetic/North']);
    tracks.(sat) = tmp;
    clear('tmp')
end

%% replicate
plot_suffix = strrep(strip(char(cfg.plot_suffix),'_'),' ','');
if ~isempty(plot_suffix)
    plot_suffix = ['_',plot_suffix];
end
direc_out = fullfile(direc,'output');
if ~exist(direc_out,'dir')
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
    flow_bg = [nan,nan];
end
try
    swap_primary = cfg.swap_primary;
catch
    swap_primary = false;
end

[phi,resnorm,E2_bg,E3_bg,v2_int,v3_int,weight0,bound] ...
    = aurogem.tools.replicate(tracks,image,xg ...
    ,flow_smoothing_window = cfg.flow_smoothing_window ...
    ,boundary_smoothing_window = cfg.boundary_smoothing_window ...
    ,show_plots = true ...
    ,save_plots = true ...
    ,direc = direc_out ...
    ,suffix = plot_suffix ...
    ,add_phi_background = cfg.add_phi_background ...
    ,fit_harmonic = cfg.fit_harmonic_function ...
    ,num_replications = cfg.replication_number ...
    ,arc_definition = cfg.arc_definition ...
    ,edge_method = "contour" ...
    ,do_rotate = true ...
    ,do_scale = true ...
    ,contour_values = cfg.contour_values ...
    ,harmonic_mask = cfg.harmonic_mask ...
    ,weighting_scale_length = wsl ...
    ,flow_bg = flow_bg ...
    ,swap_primary = swap_primary ...
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
    answer = input('Save potential map? (Y/n) ','s');
    if strcmp(answer,'Y')
        answered = true;
        happy = true;
    elseif strcmp(answer,'n')
        answered = true;
        happy = false;
    end
end

if happy
    file_phi = fullfile(direc,'potential.h5');
    fprintf('Saving %s\n',file_phi)
    h5make(file_phi,'/Driver/Potential',phi,'Convection Electric Potential' ...
        ,units='Volts',size='lxp x lyp')
    h5make(file_phi,'/Driver/ResNorm',resnorm,'Relative residual norm' ...
        ,size='1 x 2 (East, North)',note='norm(residual)/norm(interpolated)')
    h5make(file_phi,'/BackgroundE/East',E2_bg ...
        ,'Magnetic eastward background electric field',units='Volts/meter')
    h5make(file_phi,'/BackgroundE/North',E3_bg ...
        ,'Magnetic northward background electric field',units='Volts/meter')
    h5make(file_phi,'/InterpolatedFlow/East',v2_int,'Magnetic eastward interpolated flow' ...
        ,units='Meters/second',size='lxp x lyp')
    h5make(file_phi,'/InterpolatedFlow/North',v3_int,'Magnetic northward interpolated flow' ...
        ,units='Meters/second',size='lxp x lyp')
    for sat = sats
        h5make(file_phi,['/Weights/',sat],weight.(sat),['Track ',sat,' weight map'] ...
            ,size='lxp x lyp')
    end
    h5make(file_phi,'/Boundary/Primary',bound_prim,'Primary arc boundary'...
        ,units='Meters',size='lxp x 2 (east, north)')
    h5make(file_phi,'/Boundary/Secondary',bound_scnd,'Secondary arc boundary' ...
        ,units='Meters',size='lxp x 2 (east, north)')
end

close all
fclose all;
end