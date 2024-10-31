function generate_invert_table(date, glat, glon, opts)
arguments
    date (1, 1) datetime
    glat (1, 1) double {mustBeNonempty}
    glon (1, 1) double {mustBeNonempty}
    opts.version (1, 1) int32 {mustBePositive} = 3
    opts.out_direc (1, :) char {mustBeFolder} = '.'
    opts.do_acc (1, 1) logical = false
    opts.thermal_energy_ev (1, 1) double = 0;   
end

assert(not(ispc), 'Please run this function in a Linux environment.')

ver = sprintf('v%i', opts.version);
in_filename = jules.glow.generate_input(date, glat, glon, ...
    do_acc=opts.do_acc, thermal_energy_ev=opts.thermal_energy_ev);
in_path = fullfile('input', in_filename);

for airglow = [true, false]

    glow_out_path = fullfile('glow', 'output', ver);
    tmp = strsplit(in_filename, '.');
    if strcmp(opts.out_direc, '.')
        out_direc = fullfile(glow_out_path, tmp{3});
    else
        out_direc = fullfile(opts.out_direc, tmp{3});
    end
    if airglow
        out_direc = fullfile(out_direc, 'airglow');
    end
    fprintf('\nInput filename: %s\n', in_path)
    fprintf('Glow output directory: %s\n', glow_out_path)
    fprintf('Output directory: %s\n\n', out_direc)
    
    if ~isfolder(out_direc)
        mkdir(out_direc)
    end
    copyfile(fullfile('glow', in_path), out_direc, 'f')

    if airglow
        fprintf('Generating airglow inversion value.\n')
        command = fullfile('.', sprintf('glow_invert_airglow_%s.exe < %s', ver, in_path));
    else
        fprintf('Generating inversion tables.\n')
        command = fullfile('.', sprintf('glow_invert_tables_%s.exe < %s', ver, in_path));
    end
    system(sprintf('cd glow; %s', command));
    
    bins = {dir(fullfile(glow_out_path, '*.bin')).name};
    for bin = bins
        movefile(fullfile(glow_out_path, bin{1}), out_direc, 'f')
    end
end