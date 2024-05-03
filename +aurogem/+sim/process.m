function process(direc,opts)
arguments
    direc (1,:) char {mustBeFolder}
    opts.start (1,1) double {mustBeNonempty} = -1
    opts.cad (1,1) double {mustBeNonempty} = -1
    opts.stop (1,1) double {mustBeNonempty} = -1
    opts.mlon_ref (1,1) double {mustBeNonempty} = -1
    opts.alt_ref (1,1) double {mustBePositive} = 300e3
    opts.plot (1,1) logical {mustBeNonempty} = true
    opts.video (1,1) logical {mustBeNonempty} = true
    opts.vtk (1,1) logical {mustBeNonempty} = false
end

if opts.plot
    jules.sim.plot(direc...
        ,start = opts.start...
        ,stop = opts.stop...
        ,cad = opts.cad...
        ,alt_ref = opts.alt_ref...
        ,mlon_ref = opts.mlon_ref...
        )
end

if opts.video
    plots_dirs = dir(fullfile(direc,'plots'));
    for i = 3:length(plots_dirs)
        plots_name = plots_dirs(i).name;
        jules.tools.images2video(direc,fullfile('plots',plots_name))
    end
end

if opts.vtk
    gemini3d.write.vtk(direc,'njv')
end

end