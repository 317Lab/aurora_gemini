% Description:
%   Generate plots, videos, and/or vtk files of finished simulation.
%
% Example usage:
%   aurogem.sim.process('path-to-simulation')
%
% Arguments:
%   direc                   simulation directory
%   start = -1              (option) plotting start time (-1 = cfg.dtout)
%   cad = -1                (option) plotting cadence (-1 = cfg.dtout)
%   stop = -1               (option) plotting stop time (-1 = cfg.tdur)
%   alt_ref                 (option) reference altitude
%   mlon_ref  = -1          (option) reference magnetic longitude (-1 = mean(mlon))
%   do_plot = true          (option) whether to generate plots
%   do_video = true         (option) whether to generate videos
%   do_vtk = false          (option) whether to generate vtk files
%
% Contact:
%   jules.van.irsel.gr@dartmouth.edu
%
% Revisions:
%   07/23/2024  initial implementation (jvi)
%

function process(direc,opts)
arguments
    direc (1,:) char {mustBeFolder}
    opts.start (1,1) double {mustBeNonempty} = -1
    opts.cad (1,1) double {mustBeNonempty} = -1
    opts.stop (1,1) double {mustBeNonempty} = -1
    opts.alt_ref (1,1) double {mustBePositive} = 300e3
    opts.mlon_ref (1,1) double {mustBeNonempty} = -1
    opts.do_plot (1,1) logical {mustBeNonempty} = true
    opts.do_video (1,1) logical {mustBeNonempty} = true
    opts.do_vtk (1,1) logical {mustBeNonempty} = false
end

if opts.do_plot
    aurogem.sim.plot(direc...
        ,start = opts.start...
        ,stop = opts.stop...
        ,cad = opts.cad...
        ,alt_ref = opts.alt_ref...
        ,mlon_ref = opts.mlon_ref...
        )
end

if opts.do_video
    plots_dirs = dir(fullfile(direc,'plots'));
    for i = 3:length(plots_dirs)
        plots_name = plots_dirs(i).name;
        aurogem.tools.images2video(direc,fullfile('plots',plots_name))
    end
end

if opts.do_vtk
    gemini3d.write.vtk(direc,'njv')
end

end