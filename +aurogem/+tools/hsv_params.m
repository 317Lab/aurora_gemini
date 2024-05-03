% Description:
%   Generates a set of variables that can be used in pcolor for plotting
%   hsv plots, with the hue being the angle and the saturation being the
%   magntiude of the flow.
%
% Exampe usage:
%   hsv_params(v2,v3,MLAT,MLON,ALT,300e3,mean(MLON(:)),1e3)
%
% Example usage: altitude cut
%   pcolor(squeeze(MLON(alt_ref_ind,:,:)),squeeze(MLAT(alt_ref_ind,:,:)),hsv_alt)
%   shading flat
%   title('Ion Flow')
%   xlabel('Mag. Lon.')
%   ylabel('Mag. Lat.')
%   colormap(gca,hsv_alt_map)
%   clb = colorbar;
%   colormap(clb,hsv_map_clb)
%   clb.Limits = [0,1];
%   clb.Ticks = [0,1/4,1/2,3/4,1];
%   clb.TickLabels = {'W','S','E','N','W'};
%   clb.Label.String = ['v sat. at ',num2str(hsv_sat),' m/s'];
%
% Example usage: longitude cut
%   pcolor(squeeze(MLAT(:,mlon_ref_ind,:)),squeeze(ALT(:,mlon_ref_ind,:)),hsv_mlon)
%   shading flat
%   title('Ion Flow')
%   xlabel('Mag. Lon.')
%   ylabel('Alt. [km]')
%   colormap(gca,hsv_mlon_map)
%   clb = colorbar;
%   colormap(clb,hsv_map_clb)
%   clb.Limits = [0,1];
%   clb.Ticks = [0,1/4,1/2,3/4,1];
%   clb.TickLabels = {'W','S','E','N','W'};
%   clb.Label.String = ['v sat. at ',num2str(hsv_sat),' m/s'];
%
% Arguments:
%   v2          eastward flow data
%   v3          northward flow data
%   MLAT        magnetic latitude
%   MLON        magnetic longitude
%   ALT         altitude
%   alt_ref     altitude reference
%   mlon_ref    magnetic latitude reference
%   hsv_sat     hsv saturation flow level
%
% Dependencies:
%   matlab R2020a or higher
%   gemscr
%
% Contact:
%   jules.van.irsel.gr@dartmouth.edu

function [hsv_map_clb,hsv_mlon,hsv_mlon_map,hsv_alt,hsv_alt_map] = hsv_params(v2,v3,MLAT,MLON,ALT,alt_ref,mlon_ref,hsv_sat)
arguments
    v2 (:,:,:) double {mustBeNonempty}
    v3 (:,:,:) double {mustBeNonempty}
    MLAT (:,:,:) double {mustBeNonempty}
    MLON (:,:,:) double {mustBeNonempty}
    ALT (:,:,:) double {mustBeNonempty}
    alt_ref (1,1) double {mustBeNonempty}
    mlon_ref (1,1) double {mustBeNonempty}
    hsv_sat (1,1) {mustBeNonempty}
end
cnd = isequal(size(v2),size(v3),size(MLAT),size(MLON),size(ALT));
assert(cnd,'First 5 arguments must all have equal sizes.')

[~,alt_ref_ind] = min(abs(ALT(:,1,1)-alt_ref));
[~,mlon_ref_ind] = min(abs(MLON(1,:,1)-mlon_ref));

%% creating hsv arrays
% flow angles are mapped 0-1 = WSENW
% flow magnitudes are mapped 0-1 = 0-hsv_sat
% a cyclical colormap is chosen with 4 major hues for each direction
% angles are assigned to these colors
% color saturations are scaled by flow magntiude
% color values are set to 1 for hsv_sat
% for an M by N sized plot, a M*N by 3 sized colormap is created
% assigning each cell an rgb color, the first element of which
% attributing to 0, and the last to 1
% pcolor then plots linspace(0,1,M*N) reshaped to an M by N plot
ang = (1+atan2(v3,v2)./pi)/2; % 0/1 = west, 0.25 = south, 0.5 = east, etc.
mag = min(sqrt(v2.^2+v3.^2)/hsv_sat,1); % saturated flow magnitude
hsv_map_clb = jules.tools.colorcet('C2','shift',63/256); % shifted for blue is west
hsv_tmp = rgb2hsv(hsv_map_clb(round(1+255*ang),:)); % map angles to rgb to hsv
hsv_tmp(:,2) = hsv_tmp(:,2).*mag(:); % scale saturation by flow magnitude
hsv_tmp(:,3) = (hsv_tmp(:,3)-1).*mag(:)+1; % mag=1 gives initial hsv value, mag=0 gives value of 1
hsv_map = reshape(hsv2rgb(hsv_tmp),[size(MLAT),3]); % gives each flow voxel an rgb value

hsv_mlon = reshape(linspace(0,1,numel(MLAT(:,1,:))),size(squeeze(MLAT(:,1,:))));
hsv_mlon_map = reshape(hsv_map(:,mlon_ref_ind,:,:),numel(hsv_mlon),3);
hsv_alt = reshape(linspace(0,1,numel(MLAT(1,:,:))),size(squeeze(MLAT(1,:,:))));
hsv_alt_map = reshape(hsv_map(alt_ref_ind,:,:,:),numel(hsv_alt),3);
end