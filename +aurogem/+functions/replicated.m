% Description:
%   Generates gemini input fields and precipitation files from configuration and
%   grid structures. Used with aurogem.sim.replication outputs. Option
%   for "setup_functions" in "setup" namelist.
%
% Example usage:
%   aurogem.functions.replicated(cfg, xg)
%
% Arguments:
%   cfg     gemini configuration file
%   xg      gemini grid file
%
% Contact:
%   jules.van.irsel.gr@dartmouth.edu
%
% Revisions:
%   07/23/2024  initial implementation (jvi)
%

function replicated(cfg, xg)
arguments
    cfg (1, 1) struct {mustBeNonempty}
    xg (1, 1) struct {mustBeNonempty}
end

%% unpack configuration data
ymd = cfg.ymd;
UTsec0 = cfg.UTsec0;
tdur = cfg.tdur;
dtprec = cfg.dtprec;
dtE0 = cfg.dtE0;
itprec = 0:dtprec:tdur;
itE0 = 0:dtE0:tdur;
ltE0 = length(itE0);

%% unpack grid data
x2 = xg.x2(3:end-2);
x3 = xg.x3(3:end-2);
mlon = squeeze(xg.phi(end, :, 1))*180/pi;
mlat = 90-squeeze(xg.theta(end, 1, :))*180/pi;
llon = length(mlon);
llat = length(mlat);
mlon_to_x2 = griddedInterpolant(mlon, x2);
mlat_to_x3 = griddedInterpolant(mlat, x3);
[X2, X3] = ndgrid(x2, x3);

%% create east-north-time grids
[X2_prec, X3_prec, ~] = ndgrid(x2, x3, itprec);

%% generate input data maps
direc_part = fullfile(fileparts(cfg.nml), 'ext', 'precipitation.h5');
Q_part = h5read(direc_part, '/Derived/Energy/Flux');
E0_part = h5read(direc_part, '/Derived/Energy/Characteristic');
if any(isnan(Q_part(:)))
    warning('Filling NaN(s) in %s: /Derived/Energy/Flux', direc_part)
    Q_part = fillmissing2(Q_part, 'cubic');
end
if any(isnan(E0_part(:)))
    warning('Filling NaN(s) in %s: /Derived/Energy/Characteristic', direc_part)
    E0_part = fillmissing2(E0_part, 'cubic');
end

if cfg.flagdirich == 1
    warning('NaN(s) in E0_part')
    mlon_part = h5read(direc_part, '/Coordinates/Magnetic/Longitude');
    mlat_part = h5read(direc_part, '/Coordinates/Magnetic/Latitude');
    x2_part = mlon_to_x2(mlon_part);
    x3_part = mlat_to_x3(mlat_part);
elseif cfg.flagdirich == 0
    x2_part = h5read(direc_part, '/Coordinates/Magnetic/East');
    x3_part = h5read(direc_part, '/Coordinates/Magnetic/North');
end
[X2_part, X3_part] = ndgrid(x2_part, x3_part);
fQ = griddedInterpolant(X2_part, X3_part, Q_part, 'cubic');
fE0 = griddedInterpolant(X2_part, X3_part, E0_part, 'cubic');
Q = fQ(X2_prec, X3_prec);
E0 = fE0(X2_prec, X3_prec);

%% write precip files
pg.mlon = mlon;
pg.mlat = mlat;
pg.llon = llon;
pg.llat = llat;
pg.Qit = Q;
pg.E0it = E0;
pg.times = datetime(ymd) + seconds(UTsec0+itprec);
gemini3d.write.precip(pg, cfg.prec_dir)

%% write field files
if cfg.flagdirich == 1
    direc_E0 = fullfile(fileparts(cfg.nml), 'ext', 'potential.h5');
    if ~exist(direc_E0, 'file')
        error('flagdirich = 1 but not potential.h5 found')
    end
    bc = h5read(direc_E0, '/Driver/Potential');
elseif cfg.flagdirich == 0
    direc_E0 = fullfile(fileparts(cfg.nml), 'ext', 'current.h5');
    if ~exist(direc_E0, 'file')
        error('flagdirich = 0 but not current.h5 found')
    end
    bc = h5read(direc_E0, '/Driver/Current');
else
    error('Unknown flagdirich')
end
E2_bg = h5read(direc_E0, '/BackgroundE/East');
E3_bg = h5read(direc_E0, '/BackgroundE/North');

cfg_rep = gemini3d.read.config(fullfile(fileparts(cfg.nml), 'ext'));
xg_rep = gemini3d.grid.cartesian(cfg_rep);
[X2_rep, X3_rep] = ndgrid(xg_rep.x2(3:end-2), xg_rep.x3(3:end-2));
clear('xg_rep')
fprintf('Interpolating potential onto simulation grid.')
fbc = griddedInterpolant(X2_rep, X3_rep, bc, 'spline');
bc = fbc(X2, X3);

if any(isnan(bc(:)))
    warning('Filling NaN(s) in field driver')
    bc = fillmissing2(bc, 'cubic');
end

Vmaxx1it = repmat(bc, [1, 1, ltE0]); % x2 * x3 * time
Vmaxx1it(:, :, 1) = zeros(llon, llat); % no forcing in first timestep

E.flagdirich = ones(1, ltE0)*cfg.flagdirich;
E.Exit = ones(llon, llat, ltE0)*E2_bg;
E.Eyit = ones(llon, llat, ltE0)*E3_bg;
E.Vminx1it = Vmaxx1it*cfg.flagdirich;
E.Vmaxx1it = Vmaxx1it;
E.Vminx2ist = squeeze(Vmaxx1it(1, :, :))*cfg.flagdirich;
E.Vmaxx2ist = squeeze(Vmaxx1it(end, :, :))*cfg.flagdirich;
E.Vminx3ist = squeeze(Vmaxx1it(:, 1, :))*cfg.flagdirich;
E.Vmaxx3ist = squeeze(Vmaxx1it(:, end, :))*cfg.flagdirich;
E.mlon = mlon;
E.mlat = mlat;
E.llon = llon;
E.llat = llat;
E.times = datetime(ymd) + seconds(UTsec0+itE0);
gemini3d.write.Efield(E, cfg.E0_dir);

end