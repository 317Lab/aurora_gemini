function [phi_out,harm] = find_harmonic(phi_0,xdata,Edata,xg)
arguments
    phi_0 (:,:) {mustBeNumeric}
    xdata (:,2) {mustBeNumeric}
    Edata (:,2) {mustBeNumeric}
    xg (1,1) struct {mustBeNonempty}
end
global X2 X3 dx2 dx3 phi_0_tmp order %#ok<GVMIS>

% assertions
assert(isequal(size(xdata),size(Edata)),'Sizes of xdata and Edata must match.')

% unpack grid
x2 = xg.x2(3:end-2);
x3 = xg.x3(3:end-2);
if size(x2,1)~=1; x2=x2'; end
if size(x3,1)~=1; x3=x3'; end
[X2,X3] = ndgrid(x2,x3);
dx2 = xg.dx2h;
dx3 = xg.dx3h;
% Bmag = double(abs(mean(xg.Bmag,'all')));
if all([range(dx2),range(dx3)]<1) % uniform grid?
    dx2 = mean(dx2);
    dx3 = mean(dx3);
end

% unpack optimization data
phi_0_tmp = phi_0;
[~,idata(:,1)] = min(abs(x2-xdata(:,1)),[],2);
[~,idata(:,2)] = min(abs(x3-xdata(:,2)),[],2);
Edata = double(Edata);

% find initial harmonic slopes
[E2_0,E3_0] = gradient(-phi_0_tmp,dx2,dx3);
E_0 = mean(Edata)-mean([E2_0(:);E3_0(:)]);

% regularize
reg = 1e6;
Edata = Edata*reg;
X2 = X2/reg; X3 = X3/reg;
dx2 = dx2/reg; dx3 = dx3/reg;

% optimize
order = 5;
pars_0 = zeros(1,2*order);
pars_0(1:2) = -E_0;
lb = -inf(1,2*order);
ub = -lb;
optim_opts = optimoptions('lsqcurvefit' ...
    ,'Algorithm','levenberg-marquardt' ...
    ,'FunctionTolerance',1e-30 ... % lsqcurvefit does not use opt. tol.
    ,'MaxFunctionEvaluations',3e3 ...
    ,'StepTolerance',1e-30 ...
    ,'MaxIterations',1e3 ...
    ,'Display','iter' ...
    );
fun = @(pars,idata) F(pars,idata);
[pars,~,~,exitflag] = lsqcurvefit(fun,pars_0,idata,Edata,lb,ub,optim_opts);
if exitflag <= 0
    warning('LSQCURVEFIT HAS NOT REACHED LOCAL MINIMUM')
end

harm = harmonic(pars);
% harm = tools.harmonic(X2,X3,pars,order);
phi_out = phi_0 + harm;

% functions
    function f = harmonic(p)
        rho = 10;
        f = p(1)*(X2) + p(2)*(X3) + ...
            p(3)/rho^1*(X2.^2-X3.^2) + p(4)/rho^1*(X2.*X3) + ...
            p(5)/rho^2*(X3.^3-3*X2.^2.*X3) + p(6)/rho^2*(X2.^3-3*X2.*X3.^2) + ...
            p(7)/rho^3*(X2.^4-6*X2.^2.*X3.^2+X3.^4) + ...
            p(8)/rho^3*(X2.^3.*X3-X2.*X3.^3) + ...
            p(9)/rho^4*(X2.^5+5*X2.*X3.^4-10*X2.^3.*X3.^2) + ...
            p(10)/rho^4*(X3.^5+5*X2.^4.*X3-10*X2.^2.*X3.^3);
    end
    function E = F(p,idata)
        rho = 10;
        f = p(1)*(X2) + p(2)*(X3) + ...
            p(3)/rho^1*(X2.^2-X3.^2) + p(4)/rho^1*(X2.*X3) + ...
            p(5)/rho^2*(X3.^3-3*X2.^2.*X3) + p(6)/rho^2*(X2.^3-3*X2.*X3.^2) + ...
            p(7)/rho^3*(X2.^4-6*X2.^2.*X3.^2+X3.^4) + ...
            p(8)/rho^3*(X2.^3.*X3-X2.*X3.^3) + ...
            p(9)/rho^4*(X2.^5+5*X2.*X3.^4-10*X2.^3.*X3.^2) + ...
            p(10)/rho^4*(X3.^5+5*X2.^4.*X3-10*X2.^2.*X3.^3);
%         f = tools.harmonic(X2,X3,p,order);
        phi = phi_0_tmp + f;
        E = zeros(size(idata)); % Edata and idata sizes must match
        [E2,E3] = gradient(-phi',dx2,dx3);
        for i = 1:length(idata)
            E(i,:) = [E2(idata(i,2),idata(i,1)),E3(idata(i,2),idata(i,1))];
        end
    end
end
