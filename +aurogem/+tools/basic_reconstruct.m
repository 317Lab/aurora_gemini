function phi_out = basic_reconstruct(v2,v3,xg,opts)
arguments
    v2 (:,:) {mustBeNumeric}
    v3 (:,:) {mustBeNumeric}
    xg (1,1) struct {mustBeNonempty}
    opts.usepar (1,1) logical = false
    opts.ic_iteration (1,1) logical = false
    opts.decimation (1,2) int16 {mustBePositive} = [1,1]
    opts.phi_0 (:,:) {mustBeNumeric} = []
    opts.cpu_range (1,:) int16 {mustBePositive} = [6,16]
end

% load/save data from temp file in case of crash, e.g. out-of-memory
tmp_dir = [tempdir,'basic_reconstruct.tmp'];
tmp_fn = [tempname(tmp_dir),'.mat'];
load_tmp = false;
if ~isfolder(tmp_dir)
    mkdir(tmp_dir)
elseif sum([dir(tmp_dir).bytes]) ~= 0 % tmp file found?
    tmp_fn = [dir(tmp_dir).name];
    tmp_fn = fullfile(tmp_dir,tmp_fn(4:end));
    inp = input('Load data from temp file? (y/n) ','s');
    if strcmp(inp,'y')
        load_tmp = true;
    end
end

% start timer and open parallel pool
ticA = tic;
if opts.usepar
    parpool(opts.cpu_range,'IdleTimeout',24*60)
end

% check for initial potential
load_phi_0 = false;
if not(isempty(opts.phi_0))
    load_phi_0 = true;
end

% unpack grid
lx2 = xg.lx(2);
lx3 = xg.lx(3);
[X2,X3] = ndgrid(xg.x2(3:end-2),xg.x3(3:end-2));
[DX2,DX3] = ndgrid(xg.dx2h,xg.dx3h);
% Bmag = squeeze(xg.Bmag(end,:,:));
Bmag = mean(xg.Bmag,'all');
assert(isequal(size(v2),size(v3),[lx2,lx3]), ...
    'Size of position and flow arrays must match that of grid.')
if all([range(xg.dx2h),range(xg.dx3h)]<1) % uniform grid?
    F = @F_ug;
else
    F = @F_nug;
end

% determine decimation steps
d2_min = double(opts.decimation(1));
d3 = double(opts.decimation(2));
if opts.ic_iteration && not(load_phi_0)
    d2s = 2.^(5:-1:0);
    d2s = [d2s(d2s>d2_min),d2_min];
else
    d2s = d2_min;
end
if load_tmp
    load(tmp_fn,'fphi','d2')
    d2s = d2s(d2s<d2);
elseif load_phi_0
    assert(isequal(size(opts.phi_0),[lx2,lx3]), ...
        'Size of initial potential array must match that of grid.')
    fphi = griddedInterpolant(X2,X3,opts.phi_0);
else
    fphi = griddedInterpolant(X2,X3,zeros(size(X2)));
end

% main
for d2 = d2s
    ticB = tic;
    if d2 ~= d2s(end)
        fprintf('Current eastward decimation = %i\n',d2)
    else
        fprintf('Final eastward decimation = %i\n',d2)
    end

    % pack up bag-of-vectors
    Edata = zeros([size(v3(1:d2:end,1:d3:end)),2]);
    Edata(:,:,1) = -v3(1:d2:end,1:d3:end).*Bmag; % v = ExB/B^2
    Edata(:,:,2) =  v2(1:d2:end,1:d3:end).*Bmag; % Bmag < 0
    DX = zeros(size(Edata));
    DX(:,:,1) = DX2(1:d2:end,1:d3:end)*d2;
    DX(:,:,2) = DX3(1:d2:end,1:d3:end)*d3;

    % regularize data; E [V/km] = - grad [1/km] * phi [V]
    regularization = 1e3;
    Edata = Edata*regularization; % V/km
    DX = DX/regularization; % km
    fprintf('Regularized electric field range = [ %.2f -- %.2f ]\n' ...
        ,min(Edata(:)),max(Edata(:)))
    fprintf('Regularized position gradient range = [ %.2f -- %.2f ]\n' ...
        ,min(DX(:)),max(DX(:)))

    % optimize
    phi_0 = double(fphi(X2(1:d2:end,1:d3:end),X3(1:d2:end,1:d3:end)));
    clear('fphi','phi') % free up memory
    fprintf('Current number of fitting elements %i\n',numel(phi_0))
    lb = -inf(size(phi_0));
    ub = inf(size(phi_0));
    optim_opts = optimoptions('lsqcurvefit' ...
        ,'Algorithm','levenberg-marquardt' ...
        ,'UseParallel',opts.usepar ...
        ,'FunctionTolerance',1e-5 ... % lsqcurvefit does not use opt. tol.
        ,'Display','iter' ...
        );
    [phi,resnorm,~,exitflag] = lsqcurvefit(@(phi,DX) F(phi,DX),phi_0,DX,Edata ...
        ,lb,ub,optim_opts);
    if exitflag <= 0
        warning('LSQCURVEFIT HAS NOT REACHED LOCAL MINIMUM')
    end
    fphi = griddedInterpolant( ...
        X2(1:d2:end,1:d3:end), X3(1:d2:end,1:d3:end), phi ...
        );
    save(tmp_fn,'fphi','d2') % save temp data in case of memory crash
    clear('phi_0','Edata','DX') % free up memory
    recon_time = toc(ticB);
    fprintf('Reconstruction time = %.2f seconds\n',recon_time)
end

if opts.usepar
    delete(gcp('nocreate'))
end

% resample onto 'xg' grid
if isequal(opts.decimation,[1,1])
    phi_out = phi;
else
    phi_out = fphi(X2,X3);
end

total_time = toc(ticA);
fprintf('Total reconstruction time = %.2f seconds for %i elements with a residual norm of %.2f\n', ...
    total_time,numel(phi),resnorm)
delete(tmp_fn)

% least-square curve fitting function
    function E = F_ug(phi,DX)
        E = zeros(size(DX)); % DX and Edata sizes must match
        dx2 = DX(1,1,1);
        dx3 = DX(1,1,2);
        [E2,E3] = gradient(-phi',dx2,dx3);
        E(:,:,1) = E2';
        E(:,:,2) = E3';
    end
% least-square curve fitting function - non-uniform grid
    function E = F_nug(phi,DX)
        E = zeros(size(DX));
        dx2 = DX(:,1,1)';
        dx3 = DX(1,:,2);
        [E2,E3] = gradient(-phi',dx2,dx3);
        E(:,:,1) = E2';
        E(:,:,2) = E3';
    end
end
