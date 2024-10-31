function filename = generate_input(date, glat, glon, opts)
arguments
    date (1, 1) datetime
    glat (1, 1) double {mustBeNonempty}
    glon (1, 1) double {mustBeNonempty}
    opts.do_acc (1, 1) logical = false
    opts.thermal_energy_ev (1, 1) double = 0;    
end

date.Format = 'yyDDD';
glat = mod(glat + 90, 181) - 90;
glon = mod(glon, 360);
[f107, f107p, f107a, Ap] = jules.tools.activity(date);
% Ec = 1; % not used
% Q = 1000; % not used
filename = sprintf('in.invert.%s_%i', date, second(date, 'secondofday'));
if opts.do_acc
    filename = [filename, '_acc'];
end

f = fopen(fullfile('glow', 'input', filename), 'w');
fprintf(f, '%s %5i %5.2f %6.1f %6.1f %6.1f %6.1f %6.1f %1i %7.1f\n', ...
    date, second(date, 'secondofday'), glat, glon, f107a, f107, f107p, Ap, opts.do_acc, opts.thermal_energy_ev);

fclose all;