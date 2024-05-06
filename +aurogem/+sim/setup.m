function setup(direc,hpc)
arguments
    direc (1,:) char {mustBeFolder}
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
aurogem.sim.run_ic(direc)
gemini3d.model.setup(direc,direc)
if strcmp(hpc,"hec")
    aurogem.sim.pbs(direc)
elseif strcmp(hpc,"discovery")
    aurogem.sim.slurm(direc)
elseif strcmp(hpc,"both")
    aurogem.sim.pbs(direc)
    aurogem.sim.slurm(direc)
else
    warning('HPC %s not found. No batch script made.',hpc)
end