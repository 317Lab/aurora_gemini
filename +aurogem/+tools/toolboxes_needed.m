function toolboxes_needed(varargin)
    if isempty(varargin)
        filename = dbstack(1).file;
    else
        filename = varargin;
    end
    [~,req_boxes] = matlab.codetools.requiredFilesAndProducts(filename);
    disp('Required toolboxes:')
    for req = req_boxes
        if ~strcmp(req.Name,'MATLAB')
            disp(['   - ',req.Name])
        end
    end