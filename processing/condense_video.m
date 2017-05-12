function condense_video(varargin)
% The 2x2 binning on the camera returns the video doublewide (columns are
% duplicated), so we just need to take every other column

if nargin == 0
    getfile_choice = questdlg('Directly select files, or select all files within a directory and its subdirectories?',...
                     'Select Files or Folder?',...
                     'Select Folder','Select Files','Select Files');
    switch getfile_choice
        case 'Select Folder'
            path = uigetdir;
            filename_ext = getFilenames(path,{'\w*.mj2'});
        case 'Select Files'
            [filename_ext, path] = uigetfile('*.mj2', 'please choose source files','MultiSelect','on');
            if ischar(filename_ext) % In the case of a single file...
                filename_ext = {filename_ext};
            end

            %Prepend path to files
            filename_ext = cellfun(@(x) fullfile(path,x),filename_ext,'UniformOutput',false);
    end
elseif length(varargin) == 1
    filename_ext = varargin{1};
else
    sprintf('Dont know what the hell yer passing me');
end


for i=1:length(filename_ext)

    % Open VideoReader Object
    vfile = filename_ext{i};
    [~,vname,~] = fileparts(vfile);
    vr = VideoReader(vfile);

    % Error checking, want to make sure we are only trimming redundant data
    testframe = readFrame(vr);
    F1 = testframe(:,1:2:end);
    F2 = testframe(:,2:2:end);
    dframe = F1-F2;
    if sum(dframe(:)) ~= 0
        sprintf('Skipping %s, Columns are not redundant!',vname)
        continue
    else
        %sprintf('Condensing %s',vname);
        wb = waitbar(0,sprintf('Condensing %s, File %d/%d, %.1f%%',vname,i,length(filename_ext),0));
        vr.CurrentTime = 0;

        vw = VideoWriter([vfile(1:end-4),'_condensed'],'Archival');
        vw.FrameRate = vr.FrameRate;
        vw.MJ2BitDepth = 12;
        open(vw);

        while hasFrame(vr)
            frame = readFrame(vr);
            writeVideo(vw,frame(:,1:2:end));
            pos = vr.CurrentTime/vr.Duration;
            waitbar(pos,wb,sprintf('Condensing %s, File %d/%d, %.1f%%',vname,i,length(filename_ext),pos*100))
        end

        close(vw);
        delete(wb)
    end
end

    
