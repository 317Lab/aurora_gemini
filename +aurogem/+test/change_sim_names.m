root_direc = getenv('GEMINI_SIM_ROOT');
sims = {dir(fullfile(root_direc, 'swop_2023*')).name};

for sim = sims
    old_name = sim{1};
    disp(old_name)
    data = strsplit(old_name, '_');
    YMD = data{2};
    SOD = data{3};
    VV = data{4};
    FF = data{5};
    MM = data{6};
    SS = data{7};
    
    % check time
    time = datetime(YMD, "InputFormat", "uuuuMMdd");
    time = time + seconds(str2double(SOD));
    direc = fullfile(root_direc, old_name);
    cfg = gemini3d.read.config(direc);
    dt = seconds(cfg.times(end) - time);
    assert(dt == 0, 'BAD')
    
    % find FF flows
    event_file = fullfile(getenv('AUROGEM_ROOT'), 'data', 'swop', 'event_data.txt');
    event_data = readlines(event_file);
    found = false;
    for e = event_data(2:end-1)'
        dd = strsplit(e);
        event_time = datetime(dd(2), 'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss''Z');
        if seconds(event_time - time) == 0
            SD = str2double(strsplit(dd(14), ','));
            PF = str2double(strsplit(dd(15), ','));
            if all(PF == [0, 0])
                PF = [nan, nan];
            end
            found = true;
        end
    end
    assert(found, 'BAD')
    
    % find sim FF and SS
    cfg_ext = fullfile(root_direc, old_name, 'ext', 'config.nml');
    ext_data = readlines(cfg_ext);
    found_FF = false;
    found_SS = false;
    for d = ext_data'
        if contains(d, 'flow_background')
            tmp = strsplit(d);
            FF_test_vec = str2double(strsplit(tmp(end), ','));
            found_FF = true;
        elseif contains(d, 'used_tracks')
            tmp = strsplit(d);
            SS_test = char(strrep(tmp(end), "'", ''));
            SS_test = pad(SS_test, 2, 'left', 'x');
            found_SS = true;
        end
    end
    assert(found_FF, 'BAD')
    assert(found_SS, 'BAD')

    % check FF
    if all(SD == FF_test_vec)
        FF_test = 'SD';
    elseif all(PF == FF_test_vec)
        FF_test = 'PF';
    elseif all([0, 0] == FF_test_vec)
        FF_test = 'NB';
    else
        error('BAD')
    end
    assert(strcmp(FF, FF_test), 'BAD: %s neq %s in %s', FF, FF_test, old_name)
    
    % check SS
    assert( strcmp(SS, SS_test), 'BAD: %s neq %s in %s', SS, SS_test, old_name)

    % find MM
    cfg_tmp = fullfile(root_direc, old_name, 'config.nml');
    tmp_data = readlines(cfg_tmp);
    MM_test = 'UM';
    for d = tmp_data'
        if contains(d, 'flag_fang')
            MM_test = 'AM';
            break
        end
    end

    % check MM
    assert(strcmp(MM, MM_test), 'BAD: %s neq %s in %s', MM, MM_test, old_name)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % if contains('SD', data)
    %     FF = 'SD';
    % elseif contains('nobg', data)
    %     FF = 'NB';
    % else
    %     FF = 'PF';
    % end
    % if contains('unacc', data)
    %     MM = 'UM';
    % else
    %     MM = 'AM';
    % end
    % SS = pad(data{4}, 2, 'left', 'x');
    % new_name = sprintf('swop_%s_%s_%s_%s_%s_%s', YMD, SOD, VV, FF, MM, SS);
    % fprintf('%s\t->\t%s\n', pad(old_name, 35), new_name)


    % slurm_file = fullfile(root_direc, old_name, 'slurm.script');
    % lines = readlines(slurm_file);
    % id = find(strlength(lines) > 0, 1, 'last');
    % fid = fopen(slurm_file, 'w');
    % for line = lines(1:id)'
    %     new_line = strrep(line, old_name, new_name);
    %     fprintf(fid, '%s\n', new_line);
    % end
    % fclose all;
    % 
    % pbs_file = fullfile(root_direc, old_name, 'pbs.script');
    % lines = readlines(pbs_file);
    % id = find(strlength(lines) > 0, 1, 'last');
    % fid = fopen(pbs_file, 'w');
    % for line = lines(1:id)'
    %     new_line = strrep(line, old_name, new_name);
    %     fprintf(fid, '%s\n', new_line);
    % end
    % fclose all;

    % file0 = fullfile(root_direc, old_name);
    % file1 = fullfile(root_direc, new_name);
    % movefile(file0, file1)

end