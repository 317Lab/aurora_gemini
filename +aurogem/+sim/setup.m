% Description:
%   Automatically configures, sets up, and runs initial condition simulation for
%   a simulation directory. Stored in <sim_root>/ics.
%
% Example usage:
%   aurogem.sim.setup('path-to-simulation', 'uname', 'hec')
%
% Arguments:
%   direc       simulation directory
%   username    nasa user id (for using hpc = "hec")
%   hpc         high-performance computing option ("hec", "discovery", "both")
%
% Contact:
%   jules.van.irsel.gr@dartmouth.edu
%
% Revisions:
%   07/23/2024  initial implementation (jvi)
%

function setup(direc, username, hpc)
arguments
    direc (1, :) char {mustBeFolder}
    username (1, :) char {mustBeNonempty}
    hpc (1, 1) string {mustBeMember(hpc, ["hec", "discovery", "both"])}
end

potential_fn = fullfile(direc, 'ext', 'potential.h5');
current_fn = fullfile(direc, 'ext', 'current.h5');
if not(exist(potential_fn, 'file')) && not(exist(current_fn, 'file'))
    error('Please have either potential.h5 or current.h5, not both.')
elseif exist(potential_fn, 'file') && exist(current_fn, 'file')
    error('Please run aurogem.sim.replication first.')
end

mat_root = getenv('GEMINI_MAT_ROOT');
assert(~isempty(mat_root), ...
    ['Add environment variable GEMINI_MAT_ROOT directing to contents of ' ...
    'https://github.com/gemini3d/mat_gemini'])
addpath(mat_root)

% run initial condition and setup
aurogem.sim.run_ic(direc)
gemini3d.model.setup(direc, direc)

% check if input fields on gemini grid
simgrid_fn = fullfile(direc, 'inputs', 'simgrid.h5');
simgrid_fields_fn = fullfile(direc, 'inputs', 'fields', 'simgrid.h5');
simsize_fields_fn = fullfile(direc, 'inputs', 'fields', 'simsize.h5');
phi = rad2deg(h5read(simgrid_fn, '/phi'));
theta = rad2deg(h5read(simgrid_fn, '/theta'));
mlon = squeeze(phi(1, :, 1))';
mlat = 90 - squeeze(theta(1, 1, :));
mlon_fields = h5read(simgrid_fields_fn, '/mlon');
mlat_fields = h5read(simgrid_fields_fn, '/mlat');
if max(abs([mlon_fields - mlon; mlat_fields - mlat])) < 1e-3
    fprintf('%s grid matches working grid.\n', simgrid_fields_fn)
    h5write(simsize_fields_fn, '/llon', -1)
    h5write(simsize_fields_fn, '/llat', -1)
end

% generate summary plot
clm.c = 'L17'; clm.U = 'L19'; clm.p = 'D10'; clm.j = 'D1A';
colorcet = @aurogem.tools.colorcet;
ar = [1.617, 1, 1];
cfg = gemini3d.read.config(direc);
time = cfg.times(end);
time.Format = 'uuuuMMdd';
input_fn = sprintf('%s_%012.6f.h5', time, second(time, 'secondofday'));
particles_fn = fullfile(direc, 'inputs', 'particles', input_fn);
fields_fn = fullfile(direc, 'inputs', 'fields', input_fn);
simgrid_particles_fn = fullfile(direc, 'inputs', 'particles', 'simgrid.h5');

input_mlon = h5read(simgrid_particles_fn, '/mlon');
input_mlat = h5read(simgrid_particles_fn, '/mlat');
assert(all(input_mlon == h5read(simgrid_particles_fn, '/mlon')), ...
    'particles mlon and fields mlon do not match');
assert(all(input_mlat == h5read(simgrid_particles_fn, '/mlat')), ...
    'particles mlat and fields mlon do not match');
[MLON, MLAT] = ndgrid(input_mlon, input_mlat);
Qp = h5read(particles_fn, '/Qp');
E0p = h5read(particles_fn, '/E0p');
Vmaxx1it = h5read(fields_fn, '/Vmaxx1it');

figure
tiledlayout(3,1)

nexttile
pcolor(MLON, MLAT, Qp)
shading flat
clb = colorbar;
clb.Label.String = 'Qp (mW/m^2)';
colormap(gca, colorcet(clm.U))
ylabel('mlat (°)')
xticks([]);
pbaspect(ar)

nexttile
pcolor(MLON, MLAT, E0p)
shading flat
clb = colorbar;
clb.Label.String = 'E0p (eV)';
colormap(gca, colorcet(clm.c))
ylabel('mlat (°)');
xticks([]);
pbaspect(ar)

nexttile
pcolor(MLON, MLAT, Vmaxx1it)
shading flat
clb = colorbar;
if cfg.flagdirich == 1
    clb.Label.String = 'Vmaxx1it (V)';
    colormap(gca, colorcet(clm.p))
elseif cfg.flagdirich == 0
    clb.Label.String = 'Vmaxx1it (A/m^2)';
    colormap(gca, colorcet(clm.j))
end
xlabel('mlon (°)'); ylabel('mlat (°)');
pbaspect(ar)

exportgraphics(gcf, fullfile(direc, 'inputs', 'summary.png'), 'Resolution', 600)

% write batch script
if strcmp(hpc, "hec")
    aurogem.sim.pbs(direc, username)
elseif strcmp(hpc, "discovery")
    aurogem.sim.slurm(direc)
elseif strcmp(hpc, "both")
    aurogem.sim.pbs(direc, username)
    aurogem.sim.slurm(direc)
else
    warning('HPC %s not found. No batch script made.', hpc)
end