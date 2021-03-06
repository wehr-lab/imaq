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
    vid = read(vin);
    disp(sprintf('Video loaded'))
    
    disp(sprintf('casting as double & rescaling'))
    vid = double(vid)/4096;
    
    % Define ROI
    if i == 1
        nframes = size(vid,4);
        roi_mask = roipoly(vid(:,:,1,round(nframes/2)));
    end
    
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
    for j = 1:length(uq_freqs)
        freq_inds = find(freqs == uq_freqs(j))+1; % add one because of stimparams offset
        freq_frames = find(ismember(stim_idx,freq_inds));
        vid_frames = squeeze(vid(:,:,1,freq_frames));
        %mean_vid_frames(:,:,:,j,i) = vid_frames;
        roi_mask_3d = repmat(roi_mask,1,1,length(freq_frames));
        global_freq_response = mean(squeeze(vid_frames(roi_mask_3d)));
        mean_freq_intensity(:,:,j,i) = mean(vid(:,:,1,freq_frames), 4)./global_freq_response;
    end

end
all_mean_intensity = mean(mean_freq_intensity,4);
[~,max_inds] = max(all_mean_intensity, [], 3);


%all_max_inds = mean(all_max_inds,3);

% Make and save tonotopy figure
p = figure;
colormap('jet')
imagesc(max_inds)
cb = colorbar;
%set(cb,'YTick',round(log(uq_freqs)))
saveas(p, [pathn,'vid-tonotopy.png'])

% Make and save trial-averaged video
%mean_vid_frames = mean(mean_vid_frames,5);
%ordered_vid = [];
%for a = 1:size(mean_vid_frames,4)
%    ordered_vid(:,:
    












    
    








