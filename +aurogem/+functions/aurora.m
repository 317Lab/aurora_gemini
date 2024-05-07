function aurora(cfg,xg)
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

%% create east-north-time grids
[X2prec,X3prec,ITprec] = ndgrid(x2,x3,itprec);
[X2E0,X3E0,ITE0] = ndgrid(x2,x3,itE0);

%% generate input data maps
[~,~,~,E0,Q] = input_maps(X2prec,X3prec,ITprec,cfg);
[J,E2,E3,~,~] = input_maps(X2E0,X3E0,ITE0,cfg);

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
                phi(ix2,ix3,it) = dot(E2(1:ix2,1,it),dx2(1:ix2)) + dot(E3(ix2,1:ix3,it),dx3(1:ix3));
            end
        end
    end
    Vmaxx1it = phi;
else % current driven
    Vmaxx1it = J;
end

Vmaxx1it(:,:,1) = zeros(llon,llat); % no forcing in first timestep

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
% input function
    function [J,E2,E3,E0,Q] = input_maps(x2,x3,it,cfg)
        eV2mJ = 1.60218e-16;
        vc = 2.99792e8; % m/s
        qe = 1.60218e-19; % C
        me = 5.10999e5/vc^2; % eV/m^2/s^2
        p = cfg;
        x2 = x2 - p.driftE*it;
        x3 = x3 - p.driftN*it;
        c = (p.ctr_spn/2)*tanh(2*p.ctr_slp*(x2-p.ctr_pos)/p.ctr_spn); % contour
        b = bar(x2,p.bar_pos+p.bar_vel*it,p.bar_frc,p.bar_gsl); % loading bar
        d = bar(it,p.dim_del,p.dim_frc,p.dim_tim); % dimming
        J_amp = p.K_amp/p.J_wth;
        J = J_amp*(... % A/m^2
            +              sheet(x3, c+p.J_wth/2        , p.J_wth        , p.J_gsl)...
            -(1/p.J_frc) * sheet(x3, c-p.J_frc*p.J_wth/2, p.J_frc*p.J_wth, p.J_gsl)...
            );
        E2 = zeros(size(J)); % V/m
        E3 = -J_amp*p.J_gsl/(4*p.SIGP_bkg)*(... % V/m
            +(1/p.J_frc)           * log(cosh(2*(x3-c+p.J_frc*p.J_wth)/p.J_gsl))...
            -((1+p.J_frc)/p.J_frc) * log(cosh(2*(x3-c                )/p.J_gsl))...
            +                        log(cosh(2*(x3-c-p.J_wth        )/p.J_gsl))...
            );
        E0_l = (p.E_amp_l-p.E_amp_b) * sheet(x3, c+p.E_pos_l, p.E_wth_l, p.E_gsl_l);
        E0_h = (p.E_amp_h-p.E_amp_l) * sheet(x3, c+p.E_pos_h, p.E_wth_h, p.E_gsl_h);
        if p.do_bar_E0_l; E0_l = E0_l.*b; end
        if p.do_bar_E0_h; E0_h = E0_h.*b; end
        if p.do_bar_j; J = J.*b; end
        if p.do_dim_E0_l; E0_l = E0_l.*d; end
        if p.do_dim_E0_h; E0_h = E0_h.*d; end
        if p.do_dim_j; J = J.*d; end
        E0 = p.E_amp_b + E0_l + E0_h; % eV
        n_acc = sqrt(2*me/p.E_amp_h)*p.r_acc*J_amp/qe; % m^-3 - acceleration region density chosen to give at most J_acc = max(J)
        J_acc = qe*n_acc*sqrt(E0)/sqrt(2*me); % A/m^2 - acceleration current from Kjell Rönnmark (2002) 10.1029/2002JA009294
        Q = 2*E0.*J_acc*eV2mJ/qe; % mW/m^2 - using J_acc = qe * int phi_M dE
        disp(['Acceleration region density: ',num2str(n_acc,'%e'),' m^-3'])
    end

% sheet definition
    function v = sheet(x3,pos,wdth,gsl)
        v = (tanh(2*(x3-pos+wdth/2)./gsl)-tanh(2*(x3-pos-wdth/2)./gsl))/2;
    end

% bar definition
    function v = bar(x2,pos,frac,gsl)
        v = (2-frac*(1-tanh(2*(x2-pos)/gsl)))/2;
    end
end