envs = [["GEMINI_ROOT" ...
    , "GEMINI_MAT_ROOT" ...
    , "GEMINI_SCR_ROOT" ...
    , "GEMINI_SIM_ROOT"] ...
    ;["contents of https://github.com/gemini3d/gemini3d" ...
    , "contents of https://github.com/gemini3d/mat_gemini" ...
    , "contents of https://github.com/gemini3d/mat_gemini-scripts" ...
    , "simulation directory"...
    ]];

for env = envs
    env_path = getenv(env(1));
    if not(isempty(env_path))
        addpath(fullfile(env_path))
        fprintf(' %s = %s\n',env(1),env_path)
        continue
    end
    filename = fullfile('data','init',env(1));
    fid = fopen(filename,'r');
    if fid == -1
        while true
            env_path = input(sprintf('Please enter root path to the %s:\n %s = ' ...
                , env(2), env(1)));
            env_path = [filesep,fullfile(env_path)];
            if isfolder(env_path)
                break
            end
        end
        fid = fopen(filename,'w');
        fprintf(fid,'%s',env_path);
    else
        env_path = fgetl(fid);
        fprintf(' %s = %s\n',env(1),env_path)
    end
    setenv(env(1),env_path)
    addpath(fullfile(env_path))
end

fclose all;
clear('envs','env','env_path','filename','fid','ans');