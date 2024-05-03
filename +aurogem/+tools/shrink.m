function data_shrunk = shrink(data)
arguments
    data (1,1) struct
end

if isfield(data,'x1')
    data_shrunk.lx = data.lx;
    data_shrunk.x1 = data.x1;
    data_shrunk.x2 = data.x2;
    data_shrunk.x3 = data.x3;
    data_shrunk.dx1h = data.dx1h;
    data_shrunk.dx2h = data.dx2h;
    data_shrunk.dx3h = data.dx3h;
    data_shrunk.dx1b = data.dx1b;
    data_shrunk.dx2b = data.dx2b;
    data_shrunk.dx3b = data.dx3b;
    data_shrunk.alt = data.alt; % needed for load_conductance
    data_shrunk.glat = data.glat; % needed for load_conductance
    data_shrunk.glon = data.glon; % needed for load_conductance
    data_shrunk.h1 = data.h1; % needed for load_conductance
    data_shrunk.Bmag = data.Bmag;
elseif isfield(data,'Phitop')
    data_shrunk.J1 = data.J1;
    data_shrunk.J2 = data.J2;
    data_shrunk.J3 = data.J3;
    data_shrunk.v1 = data.v1;
    data_shrunk.ne = data.ne;
    data_shrunk.Te = data.Te;
    data_shrunk.Ti = data.Ti;
else
    error('Data structure not found.')
end
clear('data')
end