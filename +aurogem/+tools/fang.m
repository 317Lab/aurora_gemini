function [qtot, z] = fang(msis_lst, Q0_mW, E0_keV, E0_char_keV)
arguments
    msis_lst
    Q0_mW (1, 1) double {mustBePositive}
    E0_keV (1, 1) double {mustBePositive}
    E0_char_keV (1, 1) double {mustBePositive}
end

k = 1.38064852e-16; % cm2 g s-2 K-1
G = 6.67408e-8; % cm3 g-1 s-2;
ME = 5.9722e27; % g;
RE = 6.371e8; % cm;
Deps = 35e-3; % keV;
erg2kev = 6.24150648e8;

Q0_keV = Q0_mW / erg2kev; %

P = [[1.24616,     1.45903,    -2.42269e-1,  5.95459e-2]; ...
    [ 2.23976,    -4.22918e-7,  1.36458e-2,  2.53332e-3]; ...
    [ 1.41754,     1.44597e-1,  1.70433e-2,  6.39717e-4]; ...
    [ 2.48775e-1, -1.50890e-1,  6.30894e-9,  1.23707e-3]; ...
    [-4.65119e-1, -1.05081e-1, -8.95701e-2,  1.22450e-2]; ...
    [ 3.86019e-1,  1.75430e-3, -7.42960e-4,  4.60881e-4]; ...
    [-6.45454e-1,  8.49555e-4, -4.28581e-2, -2.99302e-3]; ...
    [ 9.48930e-1,  1.97385e-1, -2.50660e-3, -2.06938e-3]];


if isstruct(msis_lst)
    z = msis_lst.z;
    rho = msis_lst.rho;
    T = msis_lst.T;
    n = msis_lst.n;
elseif exist(msis_lst, 'file')
    i = 1;
    nids = nan(1,20);
    read = false;
    count = false;
    cnt = -2;
    lines = readlines(msis_lst)';
    for l = lines
        if count
            cnt = cnt + 1;
        end
        if strcmp(l, '')
            if read
                count = true;
            end
            read = false;
            continue
        end
        if read
            tmp = strsplit(l);
            len = str2double(tmp(1)) + 2;
            if strcmp(tmp(2), 'Height,')
                zid = str2double(tmp(1)) + 1;
            elseif strcmp(tmp(end), 'cm-3')
                nids(i) = str2double(tmp(1)) + 1;
                i = i + 1;
            elseif strcmp(tmp(2), 'Mass_density,')
                rhoid = str2double(tmp(1)) + 1;
            elseif strcmp(tmp(2), 'Temperature_neutral,')
                Tid = str2double(tmp(1)) + 1;
            end
        end
        if strcmp(l, 'Selected output parameters:')
            read = true;
        end
    end
    nids(isnan(nids)) = [];

    i = 1;
    z = nan(1, cnt);
    n = nan(length(nids), cnt);
    rho = nan(1, cnt);
    T = nan(1, cnt);
    for l = lines(end-cnt:end-1)
        tmp = strsplit(l);
        if length(tmp) == len
            z(i) = str2double(tmp(zid)) * 1e5; % cm
            rho(i) = str2double(tmp(rhoid)); % g cm-3
            T(i) = str2double(tmp(Tid)); % K
            for j = 1:length(nids)
                n(j,i) = str2double(tmp(nids(j))); % cm-3
            end
        end
        i = i + 1;
    end
else
    error('msis_list must be msis output file or matlab struct.')
end

m = rho ./ sum(n, 1); % g
g = G * ME ./ (RE + z).^2; % cm s-2
H = k * T ./ (m .* g); % cm

lbins = 64;
qtot = 0;
% Ebin_keV_list = nan(1, lbins);
% dEbin_keV_list = nan(1, lbins);
for k = 1:lbins
    % bin_lb = -1;
    % bin_ub = 0*2 + log10(20 * E0_char_keV) + 0*log10(20 * E0_keV);
    % bin_ub = log10(20 * max(E0_char_keV, E0_keV));
    bin_lb = log10(E0_keV); % minimal energy always
    bin_ub = log10(10 * max(E0_keV, E0_char_keV));
    Ebin_keV =  (10^(bin_lb+(bin_ub-bin_lb)*(k)/(lbins-1)) + 10^(bin_lb+(bin_ub-bin_lb)*(k-1)/(lbins-1)))/2;
    dEbin_keV = (10^(bin_lb+(bin_ub-bin_lb)*(k)/(lbins-1)) - 10^(bin_lb+(bin_ub-bin_lb)*(k-1)/(lbins-1)));
    
    % Ebin_keV_list(k) = Ebin_keV;
    % dEbin_keV_list(k) = dEbin_keV;

    y = 2 / Ebin_keV * (rho .* H / 6e-6).^0.7;

    C = zeros(1, size(P, 1));
    for i = 1:size(P, 1)
        for j = 1:size(P, 2)
            C(i) = C(i) + P(i, j) * log(Ebin_keV)^(j-1);
        end
    end
    C = exp(C);

    f = C(1)*y.^C(2).*exp(-C(3)*y.^C(4)) + C(5)*y.^C(6).*exp(-C(7)*y.^C(8));

    if Ebin_keV < E0_keV
        phi_keV = 0;
    else
        phi_keV = (Q0_keV / (E0_char_keV^2 + (E0_char_keV + E0_keV)^2)) * ...
            (Ebin_keV / E0_char_keV) * exp(-(Ebin_keV - E0_keV) / E0_char_keV); % keV-1 s-1 cm-2
    end

    dQ0_keV = Ebin_keV * phi_keV * dEbin_keV; % keV s-1 cm-2
    qtot = qtot + f * dQ0_keV / Deps ./ H; % cm-3 s-1
end
% close all
% figure
% plot(log10(Ebin_keV_list), '-or'); ylabel('log Ebin (keV)')
% figure
% plot(log10(dEbin_keV_list), '-ob'); ylabel('log dEbin (keV)')