function preprocess_videos
% Basic fixes to videos:
% Column doubling
% Saturated first frame
% DFOF

% Ask for videos
[filen, pathn, ~] = uigetfile('*','MultiSelect','on');

% cd to path so we can find stimparams
cd(pathn);

% Iterate through videos
for i=1:length(filen)
    
    % Load Video
    tic;
    disp(sprintf('Loading %s, %d of %d ...',[pathn,filen{i}],i,length(filen)))
    vr = VideoReader([pathn,filen{i}]);
    vid = read(vr);
    disp(sprintf('Loaded in %.2f',toc))
    
    % Load stimparams and parse
    disp(sprintf('Loading stimparams'))
    this_filen = filen{i};
    this_filen = this_filen(1:end-4);
    params_fn = dir([this_filen,'*.mat']);
    params = load(params_fn.name);
    params = params.userdata;
    stim_idx = params.stim_frame_idx;
    events = params.camera_events;
    stimuli = params.stimuli;

    % Fix column doubling if relevant
    if size(vid,1) ~= size(vid,2)
        disp('Video not square, Fixing column doubling')
        vid = vid(:,1:2:end,:,:);
    end
    
    % Fix saturation on first frame
    vid(:,:,:,1) = mean(vid(:,:,:,:),4);
    
    % Do dfof
    disp(sprintf('Casting as double'))
    vid = double(vid);
    
    disp(sprintf('Doing Dfof'))
    resting_frames = stim_idx==0; 
    % take mean and divide
    tic
    mean_pix = mean(vid(:,:,:,resting_frames),4);
    for m=1:size(vid,4)
        vid(:,:,:,m) = vid(:,:,:,m)./mean_pix;
    end
    disp(sprintf('dfof completed in in %.2f',toc))
    
    disp(sprintf('Saving as mj2'))
    % Scale video so .5 is 1, and proportionality is preserved
    vid = vid-1;
    if max(vid(:))>=abs(min(vid(:)))
        vid = (vid/(max(vid(:))*2))+0.5;
    else
        vid = (vid/(min(vid(:))*2))+0.5;
    end
    
    vid_int = uint16(vid*4096);
    
    tic
    vw = VideoWriter([this_filen,'_preproc'],'Archival');
    vw.FrameRate = vr.FrameRate;
    vw.MJ2BitDepth = 12;
    open(vw);
    writeVideo(vw,vid_int);
    close(vw);
    disp(sprintf('mj2 write completed in in %.2f',toc))
    
    disp(sprintf('Saving as h264'))
    tic
    vw = VideoWriter([this_filen,'_preproc'],'Grayscale AVI');
    vw.FrameRate = vr.FrameRate;
    open(vw);
    writeVideo(vw,vid);
    close(vw);
    disp(sprintf('h264 write completed in in %.2f',toc))

end