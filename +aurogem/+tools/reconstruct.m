function [phi,P] = reconstruct(pos,flow,boundary,xg,options)
arguments
    pos (:,2) {mustBeNumeric}
    flow (:,2) {mustBeNumeric}
    boundary (:,2) {mustBeNumeric}
    xg (1,1) struct {mustBeNonempty}
    options.numf (1,1) int32 {mustBePositive} = 32
    options.x2inc_lim (1,1) double = Inf
    options.usepar (1,1) logical = false
    options.maxiter (1,1) int16 = 400
    options.maxfuneval (1,1) int32 = 1000
    options.Bmag (1,1) {mustBeNumeric} = 5e-05
    options.res (1,1) int16 = 1
end
global F_opts %#ok<GVMIS>

% assertions
assert(isequal(size(pos),size(flow)),'Size of position and flow arrays must match.')

% unpack bag-of-vectors
xdata = pos(1:options.res:end,:);
Edata = flip(flow(1:options.res:end,:).*[-1,1]*options.Bmag,2); % v = ExB/B^2
disp(['Size of bag = ',num2str(length(xdata))])

% unpack grid
[X2,X3] = ndgrid(double(xg.x2(3:end-2)),double(xg.x3(3:end-2)));

% interpolate boundary
Nm = options.numf; % number of basis functions
Nj = 1;
Q_0 = [zeros(Nj,1), zeros(Nj,1), zeros(Nj,1), ones(Nj,1)];
Q = lsqcurvefit(@(Q,x) Fb(Q,x),Q_0,boundary(:,1),boundary(:,2));
F_opts.Q = Q;
F_opts.Nm = Nm;
% F_opts.b = griddedInterpolant(boundary(:,1),boundary(:,2));
% F_opts.db = griddedInterpolant(boundary(1:end-1,1),diff(boundary(:,2))./diff(boundary(:,1)));

% optimize
tic
if options.usepar; parpool([6,16]); end
P_01 = max(X3(:))/1e5;
P_0 = [linspace(-P_01,P_01,Nm); ones(1,Nm); zeros(1,Nm); ones(1,Nm)]; % initial parameter matrix: [x3pos (100 km), x3sig (100 km), x2inc (kV / 1e5 km), x2amp (kV)]
lb = repmat(-[Inf Inf options.x2inc_lim Inf]',[1 Nm]); % width structure limits, TBD: somehow largest s/c sep. %%% CHANGE
ub = repmat( [Inf Inf options.x2inc_lim Inf]',[1 Nm]);
lsq_opts = optimoptions('lsqcurvefit' ...
    ,'Algorithm','levenberg-marquardt' ...
    ,'UseParallel',options.usepar ...
    ,'StepTolerance',1e-4 ...
    ,'MaxIterations',options.maxiter ...
    ,'MaxFunctionEvaluations',2*numel(P_0)*options.maxfuneval ...
    ,'Display','iter' ...
    );
[P,~,~,exitflag] = lsqcurvefit(@(P,xdata) F(P,xdata),P_0,xdata,Edata,lb,ub,lsq_opts); % optimized parameter matrix
if options.usepar; delete(gcp('nocreate')); end
recon_time = toc;

disp(['Reconstruction time = ',num2str(recon_time),' seconds'])
if exitflag <= 0
    warning('LSQCURVEFIT HAS NOT REACHED LOCAL MINIMUM')
end

phi = potential(P,X2,X3);
phi = phi - mean(phi,'all');

% boundary function
    function b = Fb(Q,x)
        Q = Q*1e5;
        b = zeros(size(x));
        for i = 1:size(Q,1)
            b = b + Q(i,1) + Q(i,2)*tanh((x-Q(i,3))/Q(i,4));
        end
    end

% boundary derivative function
    function db = Fdb(Q,x)
        Q = Q*1e5;
        db = zeros(size(x));
        for i = 1:size(Q,1)
            db = db + Q(i,2)*sech((x-Q(i,3))/Q(i,4)).^2/Q(i,4);
        end
    end

% phi basis function
    function phi = potential(P,x2,x3)
        x3pos = P(1,:)*1e5;
        x3sig = P(2,:)*1e5;
        x2inc = P(3,:)*1e-4;
        x2int = P(4,:)*1e3;
        phi = 0;
        for m = 1:F_opts.Nm
            %             b = F_opts.b(x2);
            b = Fb(F_opts.Q,x2);
            phi = phi + (x2inc(m).*x2 + x2int(m)).*exp(-((x3-x3pos(m)-b)./x3sig(m)).^2);
        end
    end

% least-square curve fitting function
    function E = F(P,xdata)
        Ni = size(xdata,1); % number of vectors in bag
        x3pos = P(1,:)*1e5; % latitudinal positions [m] (P entries are near unity)
        x3sig = P(2,:)*1e5; % latitudinal widths [m]
        x2inc = P(3,:)*1e-4; % longitudenal slope of potential ridge [V/m]
        x2amp = P(4,:)*1e3; % central amplitude of potential ridge [V]
        E = zeros(Ni,2); % electric field
        for i = 1:Ni % iterate through bag of vectors
            x2 = xdata(i,1); % current east position
            x3 = xdata(i,2); % current north position
            %             b = F_opts.b(x2);
            %             db = F_opts.db(x2);
            b = Fb(F_opts.Q,x2);
            db = Fdb(F_opts.Q,x2);
            % caluclate the elements of -grad(phi) = -sum_m grad(phi_m) (see documentation for details)
            for m = 1:F_opts.Nm % iterate through number of basis functions
                expf = exp(-((x3-x3pos(m)-b)/x3sig(m))^2);
                E(i,1) = E(i,1) + (-x2inc(m)-(2/x3sig(m)^2)*(x2inc(m)*x2+x2amp(m))*(x3-x3pos(m)-b)*db)*expf;
                E(i,2) = E(i,2) + (2/x3sig(m)^2)*(x2inc(m)*x2+x2amp(m))*(x3-x3pos(m)-b)*expf;
            end
        end
    end
end
