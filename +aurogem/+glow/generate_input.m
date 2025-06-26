function filename = generate_input(direc, date, glat, glon, opts)
arguments
    direc (1, :) char {mustBeFolder}
    date (1, 1) datetime
    glat (1, 1) double {mustBeNonempty}
    glon (1, 1) double {mustBeNonempty}
    opts.do_acc (1, 1) logical = false
    opts.thermal_energy_ev (1, 1) double = 0;    
end

date.Format = 'yyDDD';
glat = mod(glat + 90, 181) - 90;
glon = mod(glon, 360);
[f107, f107p, f107a, Ap] = aurogem.tools.activity(date);
% Ec = 1; % not used
% Q = 1000; % not used
filename = sprintf('in.invert.%s_%i', date, second(date, 'secondofday'));
if opts.do_acc
    filename = [filename, '_acc'];
end

input_direc = fullfile(direc, 'input');
if ~exist(input_direc, 'dir')
    mkdir(input_direc)
end
f = fopen(fullfile(input_direc, filename), 'w');
fprintf('Writing %s\n', fullfile(direc, 'input', filename))
fprintf(f, '%s %5i %5.2f %6.1f %6.1f %6.1f %6.1f %6.1f %1i %7.1f\n', ...
    date, second(date, 'secondofday'), glat, glon, f107a, f107, f107p, Ap, opts.do_acc, opts.thermal_energy_ev);

fclose all;