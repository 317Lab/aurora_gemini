function h5make(file,name,var,title,opts)
arguments
    file (1,:) char
    name (1,:) char
    var (:,:)
    title (1,:) char
    opts.type (1,:) char = 'double'
    opts.units (1,:) char = ''
    opts.size (1,:) char = ''
    opts.foot_alt (1,:) char = ''
    opts.note (1,:) char = ''
end
h5create(file,name,size(var),'Datatype',opts.type)
h5write(file,name,var)
h5writeatt(file,name,'Title',title)
if ~isempty(opts.units)
    h5writeatt(file,name,'Units',opts.units)
end
if ~isempty(opts.size)
    h5writeatt(file,name,'Size',opts.size)
end
if ~isempty(opts.foot_alt)
    h5writeatt(file,name,'FootpointAltitude',opts.foot_alt)
end
if ~isempty(opts.note)
    h5writeatt(file,name,'Note',opts.note)
end
end