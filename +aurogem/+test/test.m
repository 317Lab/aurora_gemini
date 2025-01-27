msis_lst = '../Jules/thesis/data/msis_fang_2010.lst';
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
msis.z = nan(1, cnt);
msis.n = nan(length(nids), cnt);
msis.rho = nan(1, cnt);
msis.T = nan(1, cnt);
for l = lines(end-cnt:end-1)
    tmp = strsplit(l);
    if length(tmp) == len
        msis.z(i) = str2double(tmp(zid)) * 1e5; % cm
        msis.rho(i) = str2double(tmp(rhoid)); % g cm-3
        msis.T(i) = str2double(tmp(Tid)); % K
        for j = 1:length(nids)
            msis.n(j,i) = str2double(tmp(nids(j))); % cm-3
        end
    end
    i = i + 1;
end

fn = '//dartfs-hpc/rc/lab/L/LynchK/public_html/Gemini3D/swop_20230210_35487_A_05/inputs/particles/20230210_35487.000000.h5';
fn_grid = '//dartfs-hpc/rc/lab/L/LynchK/public_html/Gemini3D/swop_20230210_35487_A_05/inputs/particles/simgrid.h5';
Qp = h5read(fn, '/Qp');
E0p = h5read(fn, '/E0p');
mlat = h5read(fn_grid, '/mlat');
mlon = h5read(fn_grid, '/mlon');
[MLON, MLAT] = ndgrid(mlon, mlat);
E0_char_keV = 0.490;

dec = 1;
MLON = MLON(1:dec:end, 1:dec:end);
MLAT = MLAT(1:dec:end, 1:dec:end);
Qp = Qp(1:dec:end, 1:dec:end);
E0p = E0p(1:dec:end, 1:dec:end);

rate = nan(size(Qp, 1), size(Qp, 2), 201);
for i = 1:size(Qp, 1)
    disp(i)
    for j = 1:size(Qp, 2)
        Q0_mW = Qp(i, j);
        E0_keV = E0p(i, j) / 1e3;
        [rate(i, j, :), alt] = aurogem.tools.fang(msis, Q0_mW, E0_keV, E0_char_keV);
    end
end

%%
close all
pcolor(MLON, MLAT, sum(rate, 3));
shading flat

%%
close all
tiledlayout(5,5)

nexttile
pcolor(MLON, MLAT, E0p)
shading flat
colorbar

for i = 2:25
    id = 1 + round((i-1)*(length(alt)-1)/24);
    nexttile
    pcolor(MLON, MLAT, rate(:,:,id));
    shading flat
    xlabel(sprintf('%i km', alt(id)/1e5))
    colorbar
end

%%
[x, y] = ndgrid(MLAT(1,:), alt/1e5);

close all
z = squeeze(rate(round(size(rate, 1)/2), :, :));
pcolor(x, y, z)
shading flat
colorbar
