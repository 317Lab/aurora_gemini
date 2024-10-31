function table = read(direc)
arguments
    direc (1,:) char {mustBeFolder}
end

fn_red = fullfile(direc, dir(fullfile(direc, 'I6300*.bin')).name);
fn_gre = fullfile(direc, dir(fullfile(direc, 'I5577*.bin')).name);
fn_blu = fullfile(direc, dir(fullfile(direc, 'I4278*.bin')).name);
fn_ped = fullfile(direc, dir(fullfile(direc, 'ped3d*.bin')).name);
fn_hal = fullfile(direc, dir(fullfile(direc, 'hall3d*.bin')).name);
fn_den = fullfile(direc, dir(fullfile(direc, 'edens*.bin')).name);

[~, table.params, table.nq, table.nec, table.nalt, table.qvec, ...
    table.ecvec, table.altvec] = read3d(fn_ped);

table.red = read2d(fn_red);
table.green = read2d(fn_gre);
table.blue = read2d(fn_blu);
table.sigP = read3d(fn_ped);
table.sigH = read3d(fn_hal);
table.ne = read3d(fn_den) * 1e6;

table.units = dictionary([ ...
    "qvec", ...
    "ecvec", ...
    "altvec", ...
    "red", ...
    "green", ...
    "blue", ...
    "sigP", ...
    "sigH", ...
    "ne"
    ], [ ...
    "milliwatts meter^-2", ...
    "electronvolts", ...
    "kilometers", ...
    "Rayleighs", ...
    "Rayleighs", ...
    "Rayleighs", ...
    "Siemens", ...
    "Siemens", ...
    "meters^-3", ...
    ]);

    function [data, params, nq, nec, qvec, ecvec] = read2d(filename)
        fid = fopen(filename, 'rb');
        if fid == -1
            error('Failed to open the file: %s', filename);
        end
        raw = fread(fid, 'float32');
        fclose(fid);

        params = raw(1 : 20);
        nq = raw(21);
        nec = raw(22);
        qvec = raw(23 : 23 + nq - 1);
        ecvec = raw(23 + nq : 23 + nq + nec - 1);
        data = reshape(raw(23 + nq + nec : end), [nq, nec]);
    end

    function [data, params, nq, nec, nalt, qvec, ecvec, altvec] = read3d(filename)
        fid = fopen(filename, 'rb');
        if fid == -1
            error('Failed to open the file: %s', filename);
        end
        raw = fread(fid, 'float32');
        fclose(fid);

        params = raw(1 : 20);
        nq = raw(21);
        nec = raw(22);
        nalt = raw(23);
        qvec = raw(24 : 24 + nq - 1);
        ecvec = raw(24 + nq : 24 + nq + nec - 1);
        altvec = raw(24 + nq + nec : 24 + nq + nec + nalt - 1) * 1e-5;
        data = reshape(raw(24 + nq + nec + nalt : end), [nq, nec, nalt]);
    end
end