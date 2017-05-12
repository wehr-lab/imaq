function tonotopy
% Make an estimate of tonotopy from a calcium imaging video
%
% Asks user to select a video with a corresponding -stimparams.mat in the
% same folder, generates a dfof if dfof not in the filename, then estimates the best
% frequency for each pixel - hue for frequency, value for quality of tuning
%
% Saves image in same folder, then shows it
%%%%%%%%%%%%%%%%%%%%%%%%%

% Ask for video
[filen, pathn, ~] = uigetfile('*');

% Load video 
disp(sprintf('Reading %s ...',[pathn,filen]))
vin = VideoReader([pathn,filen]);
%vid = double(read(vin));
vid = read(vin);
vid = double(vid);

% Issue where first frame is max value, setting to mean
vid(:,:,:,1) = mean(vid(:));


%vid = (vid-min(vid(:)))/(max(vid(:))-min(vid(:)));
%vid = vid-mean(vid(:));
%vid = vid/std(vid(:));

disp(sprintf('Video loaded',[pathn,filen]))

% Load stimparams and parse
params = load([pathn, filen(1:end-4), '-stimparams.mat']);
params = params.userdata;
stim_idx = params.stim_frame_idx;
events = params.camera_events;
stimuli = params.stimuli;

% TODO: Detect if missing frames and fix stim_idx
if length(stim_idx) ~= length(events)-1
    disp('Warning! Missing some frames!')
end


% Check if square, condense if not
if size(vid,1) ~= size(vid,2)
    disp('Video not square, assuming this is the column doubling problem and condensing')
    vid = vid(:,1:2:end,:,:);
end

% Check if dfof, dfof if not
if ~findstr(filen,'dfof')
    disp('Input video is raw video, doing dfof')
    
    % get index of no-stim frames
    resting_frames = stim_idx==0;
    
    % take mean and divide
    mean_pix = mean(vid(:,:,:,resting_frames),4);
    for i=1:size(vid,4)
        vid(:,:,:,i) = vid(:,:,:,i)./mean_pix;
    end
end

% Rescale again
% vid = (vid-min(vid(:)))/(max(vid(:))-min(vid(:)));

% Find unique frequencies
freqs = [];
for i = 2:length(stimuli)
    try % Some stim may be white noise, ignore them
        freqs = [freqs, stimuli(i).param.frequency];
    end
end
uq_freqs = unique(freqs);

% Iterate over unique frequencies, average dfof intensities for each
mean_freq_intensity = [];
for i = 1:length(uq_freqs)
    freq_inds = find(freqs == uq_freqs(i))+1; % add one because of stimparams offset
    freq_frames = find(ismember(stim_idx,freq_inds));
    mean_freq_intensity(:,:,i) = mean(vid(:,:,1,freq_frames), 4);
end

% Get index of max freq
[~,max_inds] = max(mean_freq_intensity, [], 3);
max_freqs = uq_freqs(max_inds);

p = figure;
colormap('jet')
imagesc(max_freqs)
cb = colorbar;
set(cb,'YTick',round(uq_freqs))
saveas(p, [pathn,filen(1:end-4),'-tonotopy.png'])











    
    








