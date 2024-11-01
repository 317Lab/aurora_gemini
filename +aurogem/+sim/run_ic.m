% Description:
%   Automatically configures, sets up, and runs initial condition simulation for
%   a simulation directory. Stored in <sim_root>/ics.
%
% Example usage:
%   aurogem.sim.run_ic('path-to-simulation')
%
% Arguments:
%   direc               simulation directory
%   lxp = 24            number of cells in x2
%   lyp = 24            number of cells in x3
%   do_setup = true     whether to run gemini3d.model.setup
%   do_run = true       whether to run simulation
%   np = 16             number of processors for simulation
%
% Contact:
%   jules.van.irsel.gr@dartmouth.edu
%
% Revisions:
%   07/23/2024  initial implementation (jvi)
%

function run_ic(direc, opts)
arguments
    direc (1, :) char {mustBeFolder}
    opts.lxp (1, 1) int32 {mustBePositive} = 24
    opts.lyp (1, 1) int32 {mustBePositive} = 24
    opts.do_setup (1, 1) logical = true
    opts.do_run (1, 1) logical = true
    opts.np (1, 1) int32 {mustBePositive} = 16
end

%% init
gem_root = getenv('GEMINI_ROOT');
mat_root = getenv('GEMINI_MAT_ROOT');
sim_root = getenv('GEMINI_SIM_ROOT');

assert(~isempty(gem_root), ...
    'Add environment variable GEMINI_ROOT directing to gemini.bin')
assert(~isempty(mat_root), ...
    ['Add environment variable GEMINI_MAT_ROOT directing to contents of ' ...
    'https://github.com/gemini3d/mat_gemini'])
assert(~isempty(sim_root), ...
    'Add environment variable GEMINI_SIM_ROOT directing to gemini simulations')

% setup for gemini matlab tools
addpath(mat_root)
setup

% check if ic exists
cfg = gemini3d.read.config(direc);
[~, fn_ic] = fileparts(cfg.eq_dir);
direc_ic = fullfile(sim_root, 'ics', fn_ic);
if exist(direc_ic, 'dir')
    fprintf('Initial condition exists: %s\n', direc_ic)
    return
else
    mkdir(direc_ic)
end

% write ic config
fid = fopen(cfg.nml, 'r');
fid_ic = fopen(fullfile(direc_ic, 'config.nml'), 'w');
do_write = false;
is_first = true;
num_nml_written = 0;
while ~feof(fid)
    line = fgetl(fid);
    if ismember(line, {'&base', '&setup', '&flags', '&files'})
        if is_first
            fprintf(fid_ic, '%s\n', line);
            is_first = false;
        else
            fprintf(fid_ic, '\n%s\n', line);
        end
        do_write = true;
        num_nml_written = num_nml_written + 1;
        continue
    elseif strcmp(line, '/')
        fprintf(fid_ic, '/\n');
        do_write = false;
    end
    if do_write
        if startsWith(line, 'ymd')
            ymd = str2double(split(regexprep(line, '[^0-9,]', ''), ', '))';
            date = datetime(ymd) - days(1);
            ymd_new = ['ymd = ', char(date, 'uuuu, M, d')];
            fprintf(fid_ic, '%s\n', ymd_new);
        elseif startsWith(line, 'tdur')
            fprintf(fid_ic, 'tdur = %i\n', 3600*24);
        elseif startsWith(line, 'dtout')
            fprintf(fid_ic, 'dtout = %i\n', 3600*4);
        elseif startsWith(line, 'alt_scale')
            fprintf(fid_ic, 'alt_scale = 6e3, 4e3, 250e3, 3e4\n');
        elseif startsWith(line, 'x2parms')
            continue
        elseif startsWith(line, 'x3parms')
            continue
        elseif startsWith(line, 'lxp')
            line = regexprep(line, '\d', '');
            fprintf(fid_ic, '%s%i\n', line, opts.lxp);
        elseif startsWith(line, 'lyp')
            line = regexprep(line, '\d', '');
            fprintf(fid_ic, '%s%i\n', line, opts.lyp);
        elseif startsWith(line, 'eq_dir')
            fprintf(fid_ic, 'nmf = 5e11\n');
        elseif startsWith(line, 'setup_functions')
            fprintf(fid_ic, 'nme = 2e11\n');
        else
            fprintf(fid_ic, '%s\n', line);
        end
    elseif num_nml_written == 4
        break
    end
end
fclose(fid);
fclose(fid_ic);

if opts.do_setup
    gemini3d.model.setup(direc_ic, direc_ic)
end

if opts.do_run
    gemini_bin = fullfile(gem_root, 'build', 'gemini.bin');
    command = sprintf('mpiexec -np %i %s %s', opts.np, gemini_bin, direc_ic);
    system(command, '-echo')
end

end