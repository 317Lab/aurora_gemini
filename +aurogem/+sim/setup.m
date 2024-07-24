% Description:
%   Automatically configures, sets up, and runs initial condition simulation for
%   a simulation directory. Stored in <sim_root>/ics.
%
% Example usage:
%   aurogem.sim.setup('path-to-simulation','uname','hec')
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

function setup(direc,username,hpc)
arguments
    direc (1,:) char {mustBeFolder}
    username (1,:) char {mustBeNonempty}
    hpc (1,1) string {mustBeMember(hpc,["hec","discovery","both"])}
end

potential_fn = fullfile(direc,'ext','potential.h5');
if ~exist(potential_fn,'file')
    error('Please run aurogem.sim.replication first.')
end

mat_root = getenv('GEMINI_MAT_ROOT');
assert(~isempty(mat_root), ...
    ['Add environment variable GEMINI_MAT_ROOT directing to contents of ' ...
    'https://github.com/gemini3d/mat_gemini'])
addpath(mat_root)

% run initial condition and setup
aurogem.sim.run_ic(direc)
gemini3d.model.setup(direc,direc)

% check if input fields on gemini grid
simgrid_fn = fullfile(direc,'inputs','simgrid.h5');
simgrid_fields_fn = fullfile(direc,'inputs','fields','simgrid.h5');
simsize_fields_fn = fullfile(direc,'inputs','fields','simsize.h5');
phi = rad2deg(h5read(simgrid_fn,'/phi'));
theta = rad2deg(h5read(simgrid_fn,'/theta'));
mlon = squeeze(phi(1,:,1))';
mlat = 90 - squeeze(theta(1,1,:));
mlon_fields = h5read(simgrid_fields_fn,'/mlon');
mlat_fields = h5read(simgrid_fields_fn,'/mlat');
if max(abs([mlon_fields - mlon; mlat_fields - mlat])) < 1e-3
    fprintf('%s grid matches working grid.\n',simgrid_fields_fn)
    h5write(simsize_fields_fn,'/llon',-1)
    h5write(simsize_fields_fn,'/llat',-1)
end

% write batch script
if strcmp(hpc,"hec")
    aurogem.sim.pbs(direc,username)
elseif strcmp(hpc,"discovery")
    aurogem.sim.slurm(direc)
elseif strcmp(hpc,"both")
    aurogem.sim.pbs(direc,username)
    aurogem.sim.slurm(direc)
else
    warning('HPC %s not found. No batch script made.',hpc)
end