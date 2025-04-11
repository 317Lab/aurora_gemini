envs = [["GEMINI_ROOT" ...
    , "GEMINI_MAT_ROOT" ...
    , "GEMINI_SCR_ROOT" ...
    , "GEMINI_SIM_ROOT" ...
    , "AUROGEM_ROOT"] ...
    ;["contents of https://github.com/gemini3d/gemini3d" ...
    , "contents of https://github.com/gemini3d/mat_gemini" ...
    , "contents of https://github.com/gemini3d/mat_gemini-scripts" ...
    , "simulation directory"...
    , "contents of https://github.com/317Lab/aurora_gemini" ...
    ]];

if not(exist(fullfile('data', 'init'), 'dir'))
    mkdir(fullfile('data', 'init'))
end

for env = envs
    filename = fullfile('data', 'init', env(1));
    fid = fopen(filename, 'r');
    if fid ~= -1
        env_path = fgetl(fid);
        if ispc
            env_path = strrep(env_path, '/', '\');
        else
            env_path = strrep(env_path, '\', '/');
        end
    else
        fclose all;
        while true
            env_path = input(sprintf('Please enter root path to the %s:\n %s = ' ...
                , env(2), env(1)), 's');
            
            env_path = fullfile(env_path);
            if isfolder(env_path)
                break
            else
                warning('%s is not a directory', env_path)
            end
        end
        fid = fopen(filename, 'w');
        fprintf(fid, '%s', env_path);
        fclose(fid);
    end
    env_path = fullfile(filesep, env_path);
    fprintf(' %s = %s\n', env(1), env_path)
    setenv(env(1), env_path)
    addpath(env_path)
end

fclose all;
clear('envs', 'env', 'env_path', 'filename', 'fid', 'ans');