function mallinckrodt(cfg,xg)
arguments
    cfg (1,1) struct {mustBeNonempty}
    xg (1,1) struct {mustBeNonempty}
end

%% unpack configuration data
ymd = cfg.ymd;
UTsec0 = cfg.UTsec0;
tdur = cfg.tdur;
dtprec = cfg.dtprec;
dtE0 = cfg.dtE0;
itprec = 0:dtprec:tdur;
itE0 = 0:dtE0:tdur;
ltE0 = length(itE0);

%% unpack grid data
x2 = xg.x2(3:end-2);
x3 = xg.x3(3:end-2);
dx2 = xg.dx2f(3:end-2);
dx3 = xg.dx3f(3:end-2);
lx2 = length(x2);
lx3 = length(x3);
mlon = squeeze(xg.phi(end,:,1))*180/pi;
mlat = 90-squeeze(xg.theta(end,1,:))*180/pi;
llon = length(mlon);
llat = length(mlat);

%% create east-north-time grid
[~,X3prec,~] = ndgrid(x2,x3,itprec);
[~,X3E0,~] = ndgrid(x2,x3,itE0);

%% generate input data maps
if strcmp(cfg.func,'caseone')
    [~,~,~,E0,Q] = caseone(X3prec);
    [J,E2,E3,~,~] = caseone(X3E0);
elseif strcmp(cfg.func,'caseone_return')
    [~,~,~,E0,Q] = caseone_return(X3prec);
    [J,E2,E3,~,~] = caseone_return(X3E0);
else
    error([cfg.func,' does not match any map functions.'])
end

%% write precip files
pg.mlon = mlon;
pg.mlat = mlat;
pg.llon = llon;
pg.llat = llat;
pg.Qit = Q;
pg.E0it = E0;
pg.times = datetime(ymd) + seconds(UTsec0+itprec);
gemini3d.write.precip(pg,cfg.prec_dir)

%% write field files
if cfg.flagdirich % flow driven
    phi = zeros(size(J));
    for ix2 = 1:lx2
        for ix3 = 1:lx3
            for it = 1:ltE0
                phi(ix2,ix3,it) = -dot(E2(1:ix2,1,it),dx2(1:ix2))-dot(E3(ix2,1:ix3,it),dx3(1:ix3));
            end
        end
    end
    Vmaxx1it = phi;
else % current driven
    Vmaxx1it = J;
end

Vmaxx1it(:,:,1) = zeros(llon,llat,1); % no forcing in first timestep

E.flagdirich = ones(1,ltE0)*cfg.flagdirich;
E.Exit = ones(llon,llat,ltE0)*cfg.ap_ExBg;
E.Eyit = ones(llon,llat,ltE0)*cfg.ap_EyBg;
E.Vminx1it = zeros(llon,llat,ltE0);
E.Vmaxx1it = Vmaxx1it;
E.Vminx2ist = zeros(llat,ltE0);
E.Vmaxx2ist = zeros(llat,ltE0);
E.Vminx3ist = zeros(llon,ltE0);
E.Vmaxx3ist = zeros(llon,ltE0);
E.mlon = mlon;
E.mlat = mlat;
E.llon = llon;
E.llat = llat;
E.times = datetime(ymd) + seconds(UTsec0+itE0);
gemini3d.write.Efield(E,cfg.E0_dir);

%% functions
% map function
    function [J,E2,E3,E0,Q] = caseone(x3)
        J = (-0.25+3.1*sheet(x3,0,24e3,10e3))*1e-6;
%         J = 2.9e-6*(sheet(x3,15e3,15e3,10e3)-sheet(x3,-15e3,24e3,10e3));
        E2 = zeros(size(x3));
        E3 = (80-45*sheet(x3,15e3,19e3,8e3))*1e-3;
        E0 = 1000+7000*sheet(x3,0,15e3,9e3);
        Q = 90*sheet(x3,0,15e3,9e3);
    end
% map function with return current
    function [J,E2,E3,E0,Q] = caseone_return(x3)
        J = 3.1e-6*(sheet(x3,0,24e3,10e3)-sheet(x3,-24e3,24e3,10e3));
%         J = 2.9e-6*(sheet(x3,15e3,15e3,10e3)-sheet(x3,-15e3,24e3,10e3));
        E2 = zeros(size(x3));
        E3 = (80-45*sheet(x3,15e3,19e3,8e3))*1e-3;
        E0 = 1000+7000*sheet(x3,0,15e3,9e3);
        Q = 90*sheet(x3,0,15e3,9e3);
    end

% sheet definition
    function v = sheet(x3,pos,wdth,gsl)
        v = (tanh(2*(x3-pos+wdth/2)./gsl)-tanh(2*(x3-pos-wdth/2)./gsl))/2;
    end
end