% Description:
%   Generate Portable Batch System script for use on Pleiades nodes at the NASA
%   High-End Computing (HEC) Program through the NASA Advanced Supercomputing
%   (NAS) Division at Ames Research Center.
%
% Example usage:
%   aurogem.sim.pbs('path-to-simulation','uname')
%
% Arguments:
%   direc                   simulation directory
%   username                nasa user id
%   num_cpus_per_node = 16  (option) number of cpus per node
%   num_nodes = 64          (option) number of nodes
%   max_hours = 5           (option) maximum walltime in hours
%   type = "bro"            (option) processor type ("bro", "has", "ivy", "san")
%
% Contact:
%   jules.van.irsel.gr@dartmouth.edu
%
% Revisions:
%   07/23/2024  initial implementation (jvi)
%

function pbs(direc,username,opts)
arguments
    direc (1,:) char {mustBeFolder}
    username (1,:) char {mustBeNonempty}
    opts.num_cpus_per_node (1,1) int32 {mustBePositive} = 16
    opts.num_nodes (1,1) int32 {mustBePositive} = 64
    opts.max_hours (1,1) double {mustBePositive} = 5
    opts.type (1,1) string {mustBeMember(opts.type,["bro","has","ivy","san"])} = "bro"
end

%% init
gem_root = getenv('GEMINI_ROOT');
mat_root = fullfile(filesep,'u',username,'mat_gemini');
scr_root = fullfile(filesep,'u',username,'mat_gemini-scripts');
sim_root = getenv('GEMINI_SIM_ROOT');
aurogem_root = fullfile(filesep,'u',username,'aurora_gemini');

assert(~isempty(gem_root), ...
    ['Add environment variable GEMINI_ROOT directing to contents of ' ...
    'https://github.com/gemini3d/gemini3d'])
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
pfe_path = fullfile(filesep,'u',username,'sims',sim_name);

script_fn = 'pbs.script';
batch_cmd = '#PBS';
mpi_cmd = 'mpiexec';
gemini_bin = fullfile(filesep,'u',username ...
    ,'gemini3d','build','gemini.bin');
matlab_cmd = sprintf('aurogem.sim.process(''%s'')',pfe_path);
wall_h = floor(opts.max_hours);
wall_m = floor((opts.max_hours-wall_h)*60);
wall_s = round(((opts.max_hours-wall_h)*60-wall_m)*60);

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
mpis_per_node = cpus_per_node;
[n_nodes_x2,n_nodes_x3] = aurogem.tools.nearest_factors(n_nodes);
[cpus_per_node_x2,cpus_per_node_x3] ...
    = aurogem.tools.nearest_factors(cpus_per_node);
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

% write pbs script
fid = fopen(fullfile(direc,script_fn),'w');

fprintf(fid,'# Command options:\n');
fprintf(fid,'%s -N %s\n', batch_cmd, sim_name);
fprintf(fid,'%s -S /bin/bash\n', batch_cmd);
fprintf(fid,'%s -q %s\n', batch_cmd, 'normal');
fprintf(fid,'%s -l select=%i:ncpus=%i:mpiprocs=%i:model=%s\n', ...
    batch_cmd, n_nodes, cpus_per_node, mpis_per_node, opts.type);
fprintf(fid,'%s -l walltime=%i:%02d:%02d\n', batch_cmd, wall_h, wall_m, wall_s);
fprintf(fid,'%s -o %s.out\n', batch_cmd, sim_name);
fprintf(fid,'%s -e %s.err\n', batch_cmd, sim_name);
fprintf(fid,'%s -V\n', batch_cmd);

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

fprintf(fid,'\n# MPI commands:\n');
fprintf(fid,'cd /nobackup/%s\n',username);

fprintf(fid,'module use /nasa/modulefiles/testing\n');
fprintf(fid,'module load gcc/13.2\n');
fprintf(fid,'module load openmpi/4.1.6-toss4-gnu\n');

fprintf(fid,'%s %s %s\n',mpi_cmd,gemini_bin,pfe_path);

fprintf(fid,'module unuse /nasa/modulefiles/testing\n');
fprintf(fid,'module unload gcc/13.2\n');
fprintf(fid,'module unload openmpi/4.1.6-toss4-gnu\n');

fprintf(fid,'\n# Other commands:\n');
fprintf(fid,'module load matlab\n');
fprintf(fid,['matlab -batch' ...
    ' "\\\naddpath(''%s'');\\\naddpath(''%s'');\\\naddpath(''%s'');' ...
    '\\\n%s;\\\nexit\\\n"'], ...
    aurogem_root,mat_root,scr_root,matlab_cmd);

fclose(fid);
