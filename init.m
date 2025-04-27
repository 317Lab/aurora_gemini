%#ok<*UNRCH>
do_default = true; % change if you want to prompt for non-default locations

envs = ...
    [["GEMINI_ROOT" ...
    , "GEMINI_MAT_ROOT" ...
    , "GEMINI_SCR_ROOT" ...
    , "AUROGEM_ROOT" ...
    , "GEMINI_SIM_ROOT" ...
    ] ...
    ;["contents of https://github.com/317Lab/gemini3d" ...
    , "contents of https://github.com/317Lab/mat_gemini" ...
    , "contents of https://github.com/317Lab/mat_gemini-scripts" ...
    , "contents of https://github.com/317Lab/aurora_gemini" ...
    , "simulation directory"...
    ] ...
    ;["gemini3d" ...
    , "mat_gemini" ...
    , "mat_gemini-scripts" ...
    , "aurora_gemini" ...
    , "sims" ...
    ]];

if do_default
    for env = envs(:, 1:4)
        filename = fullfile('data', 'init', env(1));
        env_path = fullfile(fileparts(pwd), env(3));
        fprintf(' %s = %s\n', env(1), env_path)
        setenv(env(1), env_path)
        addpath(env_path)
    end
    query_env(envs(:, 5))
    fprintf(' EDITOR = vi\n')
    setenv('EDITOR', 'vi')
    return
end

for env = envs
    query_env(env)
end

filename = fullfile('data', 'init', 'EDITOR');
fid = fopen(filename, 'r');
if fid ~= -1
    env = fgetl(fid);
else
    fclose all;
    env = input('Please enter the text editor you want MATLAB to use:\n EDITOR = ', 's');
    fid = fopen(filename, 'w');
    fprintf(fid, '%s', env);
    fclose(fid);
end
fprintf(' EDITOR = %s\n', env)
setenv('EDITOR', env)

fclose all;
clear('do_default', 'envs', 'env', 'env_path', 'filename', 'fid', 'ans');

function query_env(env)
    if not(exist(fullfile('data', 'init'), 'dir'))
        mkdir(fullfile('data', 'init'))
    end
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