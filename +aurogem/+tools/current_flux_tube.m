function tube = current_flux_tube(xg,dat,p0,r0,r1,opts)
arguments
    xg (1,1) struct
    dat (1,1) struct
    p0 (1,3) double
    r0 (1,1) double {mustBePositive}
    r1 (1,1) double {mustBePositive}
    opts.v0 (1,3) double = [1,0,0]
    opts.v1 (1,3) double = [0,1,0]
    opts.reverse (1,1) logical = false
    opts.res (1,1) int32 = 256
    opts.do_plot (1,1) logical = false
    opts.xlims (1,2) double = [-1,1]*inf
    opts.ylims (1,2) double = [-1,1]*inf
    opts.zlims (1,2) double = [-1,1]*inf
    opts.do_reverse_test (1,1) logical = false
    opts.split_factor (1,1) double {mustBePositive} = 10
    opts.outline_axis (1,1) int8 {mustBeMember(opts.outline_axis,[2,3])} = 3;
    opts.outline_res (1,1) int32 {mustBePositive} = 1;
end

% unpack grid
scl.x = 1e-3;
x = double(xg.x2(3:end-2)*scl.x);
y = double(xg.x3(3:end-2)*scl.x);
z = double(xg.x1(3:end-2)*scl.x);
[~,lbx] = min(abs(x-max(min(x),opts.xlims(1))));
[~,lby] = min(abs(y-max(min(y),opts.ylims(1))));
[~,lbz] = min(abs(z-max(min(z),opts.zlims(1))));
[~,ubx] = min(abs(x-min(max(x),opts.xlims(2))));
[~,uby] = min(abs(y-min(max(y),opts.ylims(2))));
[~,ubz] = min(abs(z-min(max(z),opts.zlims(2))));

x = x(lbx:ubx);
y = y(lby:uby);
z = z(lbz:ubz);
dx = xg.dx2h(lbx:ubx)*scl.x;
dy = xg.dx3h(lby:uby)*scl.x;
dz = xg.dx1h(lbz:ubz)*scl.x;
[Xm,Ym,Zm] = meshgrid(x,y,z);
[dXm,dYm] = meshgrid(dx,dy);
dA = dXm.*dYm;
lz = length(z);

gap_tolerance = min([dx;dy;dz]);
split_tube_tolerance = opts.split_factor*pi*r0*r1/opts.res;

% unpack data
jx = double(permute(dat.J2(lbz:ubz,lbx:ubx,lby:uby),[3,2,1]))/scl.x^2;
jy = double(permute(dat.J3(lbz:ubz,lbx:ubx,lby:uby),[3,2,1]))/scl.x^2;
jz = double(permute(dat.J1(lbz:ubz,lbx:ubx,lby:uby),[3,2,1]))/scl.x^2;

% starting curve
v0 = opts.v0;
v1 = opts.v1;
v0 = v0/norm(v0);
v1 = v1-dot(v0,v1)*v0;
v1 = v1/norm(v1);
t = repmat(linspace(0,1,opts.res)',1,3);
c0 = p0 + r0*cos(2*pi*t).*v0 + r1*sin(2*pi*t).*v1;

% streamlines
dir = 1-2*opts.reverse;
step = 1e-2;
maxvert = 2e6;
verts = stream3(Xm,Ym,Zm,dir*jx,dir*jy,dir*jz ...
    ,c0(:,1),c0(:,2),c0(:,3),[step maxvert]);
points = vertcat(verts{:});
c1 = cell2mat(cellfun(@(v) v(end,:)',verts,'UniformOutput',0))';
c1_all = c1;

% reverse test
if opts.do_reverse_test
    for cid = 1:length(c0)
        maxvert = length(cell2mat(verts(cid)));
        vert = cell2mat(stream3(Xm,Ym,Zm,-dir*jx,-dir*jy,-dir*jz ...
            ,c1(cid,1),c1 ...
            (cid,2),c1(cid,3),[step maxvert+4]));
        delta = norm(abs(vert(end,:) - c0(cid,:)));
        if delta > gap_tolerance
            warning('Reverse test failed.')
            break
        end
    end
end

% check for split tube
c1(abs(vecnorm(diff(c1)')) > split_tube_tolerance) = nan;
if ~any(isnan(c1))
    split = false;
else
    split = true;
    for shift = 0:length(c1) % rotate until first curve in front
        if not(isnan(c1(1,1))) && isnan(c1(end,1))
            break
        end
        c1 = circshift(c1,1);
    end
    nan_ids = find(isnan(c1));
    split_ids = [1; nan_ids(diff(nan_ids)>3+1); opts.res];
    c1_tmp = cell(1,length(split_ids)-1);
    for i = 1:length(c1_tmp)
        c = c1(split_ids(i):split_ids(i+1),:);
        c(isnan(c),:) = [];
        curve_closed = [c;c(1,:)];
        if length(c) > 10
            c1_tmp(i) = mat2cell(curve_closed,length(curve_closed),3);
        else
            c1_tmp(i) = [];
        end
    end
    c1 = c1_tmp;
    verts = circshift(verts,shift);
    verts(nan_ids) = [];
end

% check flux accuracy
do_flux0 = true;
do_flux1 = true;
if range(c0(:,3)) > min(dz)
    do_flux0 = false;
    warning('Inaccurate influx due to nonhorizontal start curve.')
end
if split
    for i = 1:length(c1)
        c = cell2mat(c1(i));
        if range(c(:,3)) > min(dz)
            do_flux1 = false;
            warning('Inaccurate outflux due to nonhorizontal start curve %i.',i)
        end
    end
else
    if range(c1(:,3)) > min(dz)
        do_flux1 = false;
        warning('Inaccurate outflux due to nonhorizontal start curve.')
    end
end

% calculate fluxes
in = false(size(dXm));
out = false(size(dXm));
if do_flux0
    in = inpolygon(Xm(:,:,end),Ym(:,:,end),c0(:,1),c0(:,2));
end
if do_flux1
    if ~split
        out = inpolygon(Xm(:,:,end),Ym(:,:,end),c1(:,1),c1(:,2));
    else
        for i = 1:length(c1)
            c = cell2mat(c1(i));
            out = out | inpolygon(Xm(:,:,end),Ym(:,:,end),c(:,1),c(:,2));
        end
    end
end
flux0 = -dir*sum(squeeze(jz(:,:,end)).*dA.*in,'all');
flux1 = dir*sum(squeeze(jz(:,:,end)).*dA.*out,'all');

% generate outline and inline
ax_id = opts.outline_axis;
if ax_id == 3
    ax = z;
else
    ax = y;
end
ax = ax(1:opts.outline_res:end);
lax = length(ax);
outline = nan(2*lax-4,3);
inline = nan(2*lax-4,3);
for id = 1:lax-2
    ax0 = ax(id);
    ax1 = ax(id+1);
    layer = points(points(:,ax_id) > ax0 & points(:,ax_id) < ax1,:);
    if size(layer,1) ~= 0
        [~,min_id] = min(layer(:,5-ax_id));
        [~,max_id] = max(layer(:,5-ax_id));
        outline(lax-id-1,:) = layer(max_id,:);
        outline(lax-2+id,:) = layer(min_id,:);

        layer_sorted = sort(layer',5-ax_id)'; %#ok<UDIM>
        [gap,in_id] = max(diff(layer_sorted(:,5-ax_id)));
        if gap > gap_tolerance
            inline(lax-id-1,:) = layer_sorted(in_id,:);
            inline(lax-2+id,:) = layer_sorted(in_id+1,:);
        end
    end
end
outline = outline(not(isnan(outline(:,1))),:);
inline = inline(not(isnan(inline(:,1))),:);
outline(vecnorm(diff(outline)') > 10) = nan;

% generate tube
tube.vertices = verts;
tube.caps.start = c0;
tube.caps.end = c1;
tube.caps.end_all = c1_all;
tube.points = points;
tube.flux.in = flux0;
tube.flux.out = flux1;
tube.flux.area.in = in;
tube.flux.area.out = out;
tube.outline = outline;
tube.inline = inline;

% plotting
if opts.do_plot
    close all
    reset(0)
    set(0,'defaultLineLineWidth',2)

    fac = -squeeze(jz(:,:,end));
    fac(in | out) = nan;
    fac = repmat(fac,[1,1,lz]);

    colorcet = @jules.tools.colorcet;
    lim.x = [min(x),max(x)];
    lim.y = [min(y),max(y)];
    lim.z = [min(z),max(z)];
    lim.j = [-1,1]*quantile(abs(fac(:)),0.95);

    hold on
    title(sprintf('flux in/out ratio = %f',flux0/flux1))
    stl = streamline(verts);
    plot3(c0(:,1),c0(:,2),c0(:,3),'g')
    plot3(c0(:,1),c0(:,2),ones(size(c0(:,3)))*lim.z(1),':g')
    if ~split
        plot3(c1(:,1),c1(:,2),c1(:,3),'r')
        plot3(c1(:,1),c1(:,2),ones(size(c1(:,3)))*lim.z(1),':r')
    else
        for i=1:length(c1)
            c = cell2mat(c1(i));
            plot3(c(:,1),c(:,2),c(:,3),'r')
            plot3(c(:,1),c(:,2),ones(size(c(:,3)))*lim.z(1),':r')
        end
    end
    slice(Xm,Ym,Zm,fac,[],[],lim.z(1))
    shading flat
    colormap(colorcet('D1A'))
    xlim(lim.x); ylim(lim.y); zlim(lim.z); clim(lim.j)
    view([45,30])
    pbaspect([range(x),range(y),range(z)])
    set(stl,'Color',[0,0,0,0.1],'LineWidth',0.5)
end
end