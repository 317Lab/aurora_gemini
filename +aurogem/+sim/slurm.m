% Description:
%   Generate Slurm workload manager script.
%
% Example usage:
%   aurogem.sim.slurm('path-to-simulation')
%
% Arguments:
%   direc                       simulation directory
%   num_cpus_per_node = 64      (option) number of cpus per node
%   num_nodes = 1               (option) number of nodes
%   max_hours = 5               (option) maximum walltime in hours
%
% Contact:
%   jules.van.irsel.gr@dartmouth.edu
%
% Revisions:
%   07/23/2024  initial implementation (jvi)
%

function slurm(direc,opts)
arguments
    direc (1,:) char {mustBeFolder}
    opts.num_cpus_per_node (1,1) int32 {mustBePositive} = 64
    opts.num_nodes (1,1) int32 {mustBePositive} = 1
    opts.max_hours (1,1) double {mustBePositive} = 5
end

%% init
gem_root = getenv('GEMINI_ROOT');
mat_root = getenv('GEMINI_MAT_ROOT');
scr_root = getenv('GEMINI_SCR_ROOT');
sim_root = getenv('GEMINI_SIM_ROOT');
aurogem_root = getenv('AUROGEM_ROOT');

assert(~isempty(gem_root), ...
    ['Add environment variable GEMINI_ROOT directing to contents of ' ...
    'https://github.com/gemini3d/gemini3d'])
assert(~isempty(mat_root), ...
    ['Add environment variable GEMINI_MAT_ROOT directing to contents of ' ...
    'https://github.com/gemini3d/mat_gemini'])
assert(~isempty(scr_root), ...
    ['Add environment variable GEMINI_SCR_ROOT directing to contents of ' ...
    'https://github.com/gemini3d/mat_gemini-scripts'])
assert(~isempty(sim_root), ...
    'Add environment variable GEMINI_SIM_ROOT directing to gemini simulations')
assert(~isempty(aurogem_root), ...
    ['Add environment variable AUROGEM_ROOT directing to contents of ' ...
    'https://github.com/317Lab/aurora_gemini'])

if any(direc(end)=='/\')
    direc = direc(1:end-1);
end
[~,sim_name] = fileparts(direc);
direc = fullfile(sim_root,sim_name); % make absolute path

script_fn = 'slurm.script';
batch_cmd = '#SBATCH';
mpi_cmd = 'mpiexec';
gemini_bin = fullfile(gem_root,'build','gemini.bin');
matlab_cmd = sprintf('aurogem.sim.process(''%s'')',direc);

% check even division of cpus into horizontal cells
grid_size_fn = fullfile(direc,'inputs','simsize.h5');
lx1 = double(h5read(grid_size_fn,'/lx1'));
lx2 = double(h5read(grid_size_fn,'/lx2'));
lx3 = double(h5read(grid_size_fn,'/lx3'));
total_cells = num2str(lx1 * lx2 * lx3);
for p = 1:ceil(length(total_cells)/3)-1
    total_cells = [total_cells(1:end-4*p+1),',',total_cells(end-4*p+2:end)];
end

% check even division of cpus/nodes into horizontal cells
n_nodes = opts.num_nodes;
cpus_per_node = opts.num_cpus_per_node;
[n_nodes_x2,n_nodes_x3] = aurogem.tools.nearest_factors(n_nodes);
[cpus_per_node_x2,cpus_per_node_x3] = aurogem.tools.nearest_factors(cpus_per_node);
n_cells = lx2 * lx3;
n_cpus = double(n_nodes * cpus_per_node);
cells_per_cpu = n_cells / n_cpus;
cells_per_node_x2 = lx2 / n_nodes_x2;
cells_per_node_x3 = lx3 / n_nodes_x3;
cells_per_cpu_x2 = cells_per_node_x2 / cpus_per_node_x2;
cells_per_cpu_x3 = cells_per_node_x3 / cpus_per_node_x3;

assert(round(cells_per_cpu) == cells_per_cpu, ...
    'Number of cells (%i) do not divide evenly into number of cells (%i)', ...
    n_cells, n_cpus)
assert(round(cells_per_node_x2) == cells_per_node_x2, ...
    'Cells in x2 (%i) does not divide into number of nodes in x2 (%i)', ...
    lx2, n_nodes_x2)
assert(round(cells_per_node_x3) == cells_per_node_x3, ...
    'Cells in x3 (%i) does not divide into number of nodes in x3 (%i)', ...
    lx3, n_nodes_x3)
assert(round(cells_per_cpu_x2) == cells_per_cpu_x2, ...
    'Cells in x2 (%i) does not divide into number of cpus in x2 (%i)', ...
    cells_per_node_x2 * n_nodes_x2, cpus_per_node_x2 * n_nodes_x2)
assert(round(cells_per_cpu_x3) == cells_per_cpu_x3, ...
    'Cells in x3 (%i) does not divide into number of cpus in x3 (%i)', ...
    cells_per_node_x3 * n_nodes_x3, cpus_per_node_x3 * n_nodes_x3)

% write batch script
fid = fopen(fullfile(direc,script_fn),'w');

fprintf(fid,'#!/bin/bash\n');

fprintf(fid,'\n# Command options:\n');
fprintf(fid,'%s -J %s\n',batch_cmd,sim_name);
fprintf(fid,'%s -o %s%s%%j.log\n',batch_cmd,direc,filesep);
fprintf(fid,'%s -e %s%s%%j.err\n',batch_cmd,direc,filesep);
fprintf(fid,'%s --mail-type=END,FAIL\n',batch_cmd);
fprintf(fid,'%s --time=%i\n',batch_cmd,opts.max_hours*60);
fprintf(fid,'%s --nodes=%i\n',batch_cmd,n_nodes);
fprintf(fid,'%s --ntasks-per-node=%i\n',batch_cmd,cpus_per_node);

fprintf(fid,'\n# CPU divisions:\n');
for id = [1,fid]
    fprintf(id,'# Number of cells = %i x %i x %i = %s\n', ...
        lx1, lx2, lx3, total_cells);
    fprintf(id,'# Number of CPUs = %i\n', n_cpus);
    fprintf(id,'# Number of cells per cpu = %i\n', cells_per_cpu);
    fprintf(id,'# Number of nodes in x2 = %i\n', n_nodes_x2);
    fprintf(id,'# Number of nodes in x3 = %i\n', n_nodes_x3);
    fprintf(id,'# Number of cells per node in x2 = %i\n', ...
        cells_per_node_x2);
    fprintf(id,'# Number of cells per node in x3 = %i\n', ...
        cells_per_node_x3);
    fprintf(id,'# Number of CPUs per node in x2 = %i\n', cpus_per_node_x2);
    fprintf(id,'# Number of CPUs per node in x3 = %i\n', cpus_per_node_x3);
    fprintf(id,'# Number of cells per CPU in x2 = %i\n', cells_per_cpu_x2);
    fprintf(id,'# Number of cells per CPU in x3 = %i\n', cells_per_cpu_x3);
end

fprintf(fid,'\n# Commands to run\n');
fprintf(fid,'%s %s %s;\\\n',mpi_cmd,gemini_bin,direc);
fprintf(fid,['matlab -batch' ...
    ' "\\\naddpath(''%s'');\\\naddpath(''%s'');\\\naddpath(''%s'');' ...
    '\\\n%s;\\\nexit\\\n"'], ...
    aurogem_root,mat_root,scr_root,matlab_cmd);

fclose(fid);
end