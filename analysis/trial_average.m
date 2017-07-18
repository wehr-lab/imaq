function tonotopy_multi
% Make an estimate of tonotopy from a calcium imaging video
% Multi - average multiple shorter tuning curves b/c matlab sucks so bad.
%
% Asks user to select a video with a corresponding -stimparams.mat in the
% same folder, generates a dfof if dfof not in the filename, then estimates the best
% frequency for each pixel - hue for frequency, value for quality of tuning
%
% Saves image in same folder, then shows it
%%%%%%%%%%%%%%%%%%%%%%%%%

% Ask for videos
[filen, pathn, ~] = uigetfile('*','MultiSelect','on');
cd(pathn);

mean_freq_intensity = [];
for i=1:length(filen)

    % Load video 
    disp(sprintf('Reading %s, %d of %d ...',[pathn,filen{i}],i,length(filen)))
    vin = VideoReader([pathn,filen{i}]);
    %vid = double(read(vin));
    vid = read(vin);
    disp(sprintf('Video loaded'))
    
    disp(sprintf('casting as double & rescaling'))
    vid = double(vid)/4096;
    
    % Load stimparams and parse
    disp(sprintf('Loading stimparams'))
    this_filen = filen{i};
    this_filen = this_filen(1:3);
    params_fn = dir([this_filen,'*.mat']);
    params = load(params_fn.name);
    params = params.userdata;
    stim_idx = params.stim_frame_idx;
    events = params.camera_events;
    stimuli = params.stimuli;
    
    % Find unique frequencies
    disp('Doing mean intensity by frequency')
    freqs = [];
    for k = 2:length(stimuli)
        try % Some stim may be white noise, ignore them
            freqs = [freqs, stimuli(k).param.frequency];
        end
    end
    uq_freqs = unique(freqs);

    % Iterate over unique frequencies, find average response to each for
    % scaling
    these_frames = [];
    vid_frame_stack = [];
    for j = 1:length(uq_freqs)
        freq_inds = find(freqs == uq_freqs(j))+1; % add one because of stimparams offset
        freq_frames = find(ismember(stim_idx,freq_inds));
        vid_frames = squeeze(vid(:,:,1,freq_frames));
        for x = 1:5:length(freq_frames)
            vid_frame_stack(:,:,:,x) = vid_frames(:,:,x:x+4);
        end
        vid_frames = mean(vid_frame_stack,4);
        %vid_frames(1:20,1:20,:) = j/length(uq_freqs);
        these_frames = cat(3,these_frames,vid_frames);
    end
    all_frames(:,:,:,i) = these_frames;

end

mean_frames = mean(all_frames,4);
mean_frames = (mean_frames-min(mean_frames(:)))/(max(mean_frames(:))-min(mean_frames(:)));
z = 1;
for x=1:5:size(mean_frames,3)
    mean_frames(1:20,1:20,x:x+4) = z/length(uq_freqs);
    z = z+1;
end

% Make and save trial-averaged video
mean_frames_2(:,:,1,:) = mean_frames;
mean_frames = mean_frames_2;

vid_int = uint16(mean_frames*4096);
vw = VideoWriter('trial_average','Archival');
vw.FrameRate = vin.FrameRate;
vw.MJ2BitDepth = 12;
open(vw);
writeVideo(vw,vid_int);
close(vw);

vw = VideoWriter('trial_average','Grayscale AVI');
vw.FrameRate = vin.FrameRate;
open(vw);
writeVideo(vw,mean_frames);
close(vw);













    
    








