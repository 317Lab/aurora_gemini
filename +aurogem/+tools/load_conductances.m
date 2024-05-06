% Description:
%   Loads Pedersen and Hall conductivities/conductances from a .mat file in
%   a Gemini run subfolder if it exists. If not, it either reads them dirctly
%   from output files or reconstructs them and saves them in the that subfolder.
%
% Example usage:
%   [sigP,sigH,SIGP,SIGH] = load_conductances(direc,time,dat,cfg,xg)
%
% Arguments:
%   direc   gemini run directory
%   time    run frame datetime
%   dat     run frame data structure
%   cfg     run configuration data structure
%   xg      run grid data structure
%
% Dependencies:
%   matlab R2020a or higher
%   gemscr
%
% Contact:
%   jules.van.irsel.gr@dartmouth.edu

function [sigP,sigH,SIGP,SIGH] = load_conductances(direc,time,dat,cfg,xg)
arguments
    direc (1,:) char {mustBeFolder}
    time (1,1) datetime {mustBeNonempty}
    dat (1,1) struct {mustBeNonempty}
    cfg (1,1) struct {mustBeNonempty}
    xg (1,1) struct {mustBeNonempty}
end

h5_file = fullfile(fileparts(cfg.nml),gemini3d.datelab(time)+".h5");
try
    fprintf('Loading conductances output file... ')
    sigP = h5read(h5_file,'/sigPall');
    sigH = h5read(h5_file,'/sigHall');
    dX1 = repmat(xg.dx1h,1,xg.lx(2),xg.lx(3));
    SIGP = squeeze(sum(sigP.*dX1,1));
    SIGH = squeeze(sum(sigH.*dX1,1));
    fprintf('DONE\n')
catch
    fprintf('FAILED\n')
    folder = 'conductances';
    if ~exist(fullfile(direc,folder),'dir')
        mkdir(direc,folder);
    end
    filename_prefix = char(gemini3d.datelab(time));
    path = fullfile(direc,folder,[filename_prefix,'_sig.mat']);
    if exist(path,'file') == 2
        fprintf('Loading conductances from %s\n',path)
        load(path,'sigP','sigH','SIGP','SIGH')
    else
        fprintf('Calculating conductances (expensive!)...\n')
        [sigP,sigH,~,SIGP,SIGH,~,~] = gemscr.postprocess.conductivity_reconstruct(time,dat,cfg,xg);
        save(path,'sigP','sigH','SIGP','SIGH')
    end
end
end