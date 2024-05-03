% Combines 2 or more images in a folder to a video.
%
% Example usage:
%   images2video('..\public_html\Gemini3D\maeve_4','multipanel_auto')
%
% Arguments:
%   direc               gemini run directory
%   img_direc           subdirectory containing 2 or more images
%   ext = 'avi'         video extension ('avi','mp4','gif')
%   name = 'auto'       (option) video filename ('auto' has vid_name = img_direc)
%   img_ext = 'png'     (option) image extension
%   compressed = 1      (option) option for making an uncompressed avi
%   quality = 100       (option) 0-100 quality value if compressed
%   framerate = 5       (option) video framerate
% 
% Dependencies:
%   matlab R2020b or higher
%
% Contact: jules.van.irsel.gr@dartmouth.edu

function images2video(direc,img_direc,ext,options)
    arguments
        direc (1,:) char {mustBeFolder}
        img_direc (1,:) char {mustBeNonempty}
        ext (1,:) char {mustBeMember(ext,["avi","mp4","gif"])} = "avi"
        options.vid_name (1,:) char {mustBeNonempty} = 'auto'
        options.img_ext (1,:) char {mustBeNonempty} = 'png'
        options.compressed (1,1) logical {mustBeNonempty} = 1
        options.quality (1,1) double {mustBeNonempty} = 100
        options.framerate (1,1) double {mustBeNonempty} = 5
    end
    
    ext = char(ext);
    path = fullfile(direc,img_direc);
    img_names = dir(fullfile(path,['*.',options.img_ext]));
    if length(img_names) < 2
        error(['Fewer than 2 image files with extension ',options.img_ext,' in ',path,'.'])
    end
    
    if strcmp(options.vid_name,'auto')
        [~,options.vid_name] = fileparts(strip(img_direc,filesep));
    end
    
    if options.compressed
        if strcmp(ext,'avi')
            profile = 'Motion JPEG AVI';
        elseif strcmp(ext,'mp4')
            profile = 'MPEG-4';
        end
    else
        if ~strcmp(ext,'avi')
            warning('Changing video file extension to ".avi" for uncompressed format.')
        end
        ext = 'avi';
        profile = 'Uncompressed AVI';
        options.vid_name = [options.vid_name,'_uncompressed'];
    end
    
    if strcmp(ext,'gif')
        for i = 1:length(img_names)
            im = imread(fullfile(path,img_names(i).name));
            [imind,cm] = rgb2ind(im,256);
            if i == 1
                imwrite(imind,cm,fullfile(path,[options.vid_name,'.',ext])...
                    ,'gif','Loopcount',inf,'DelayTime',1/options.framerate);
            else
                imwrite(imind,cm,fullfile(path,[options.vid_name,'.',ext])...
                    ,'gif','WriteMode','append','DelayTime',1/options.framerate);
            end
        end
    else
        writerObj = VideoWriter(fullfile(path,[options.vid_name,'.',ext]),profile);
        writerObj.FrameRate = options.framerate;
        if options.compressed
            writerObj.Quality = options.quality;
        end
        open(writerObj)

        for i = 1:length(img_names)
            im = imread(fullfile(path,img_names(i).name));
            writeVideo(writerObj,im)
        end

        close(writerObj)
    end
end

