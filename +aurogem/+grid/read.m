function xg = read(direc, opts)
arguments
    direc (1, :) char {mustBeFolder}
    opts.all (1, 1) logical = false
end
cfg = gemini3d.read.config(direc);
xg_file = fullfile(direc, cfg.indat_grid);
xg.x1 = h5read(xg_file, '/x1');
xg.x2 = h5read(xg_file, '/x2');
xg.x3 = h5read(xg_file, '/x3');
xg.dx1h = h5read(xg_file, '/dx1h');
xg.dx2h = h5read(xg_file, '/dx2h');
xg.dx3h = h5read(xg_file, '/dx3h');
xg.Bmag = h5read(xg_file, '/Bmag');
if opts.all
   xg.glat = h5read(xg_file, '/glat'); 
   xg.glon = h5read(xg_file, '/glon'); 
   xg.alt = h5read(xg_file, '/alt'); 
   xg.phi = h5read(xg_file, '/phi'); 
   xg.theta = h5read(xg_file, '/theta'); 
end