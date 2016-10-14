function varargout =imageGUI(varargin)
% imageGUI MATLAB code for imageGUI.fig
%      IMAGEGUI, by itself, creates a new IMAGEGUI or raises the existing
%      singleton*.
%
%      H =IMAGEGUI returns the handle to a new IMAGEGUI or the handle to
%      the existing singleton*.
%
%     IMAGEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGEGUI.M with the given input arguments.
%
%     IMAGEGUI('Property','Value',...) creates a new IMAGEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imageGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imageGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imageGUI

% Last Modified by GUIDE v2.5 10-Sep-2013 13:09:40

% Begin initialization code - DO NOT EDIT
global user
global pref
if ~nargin
    [ok, user]=Login;
    if ~ok
        return
    end
end


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @imageGUI_OpeningFcn, ...
    'gui_OutputFcn',  @imageGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

if iscell(user) user=user{:};end
    h=findobj('Tag', 'User');
    set(h, 'string', user);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imageGUI_OpeningFcn(hObject, eventdata, handles, varargin)
global pref;
global user;
Prefs(user);
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and User data (see GUIDATA)
% varargin   command line arguments to imageGUI (see VARARGIN)

% Choose default command line output for imageGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
set(0,'DefaultFigureWindowStyle','normal') 

% UIWAIT makes imageGUI wait for User response (see UIRESUME)
% uiwait(handles.figure1);

Message('Did you open Sapera CamExpert and, in the Camera Link Serial Command window, and enter "set 80000" (Or however long your exposure should be given framerate)?', handles)
InitCamera(hObject, handles)
CreateDataDir(handles)

% Init sound card
OutputDeviceID          = pref.dev_id; % Sound card device #
userdata                = get(handles.figure1, 'userdata');
userdata.OutputDeviceID = OutputDeviceID;
set(handles.figure1, 'userdata', userdata);
InitSoundOut(handles)

% Load recent stimuli & init most recent
cd(fileparts(which(mfilename)));
load('RecentStimulusProtocols.mat')
stimuli=newStimList(1).stimuli;
userdata=get(handles.figure1, 'userdata');
userdata.StimList=newStimList;
set(handles.figure1, 'userdata', userdata);
InitStim(stimuli, handles);

% For live video previewing, need to add path with helper function
utilpath = fullfile(matlabroot, 'toolbox', 'imaq', 'imaqdemos', 'helper');
addpath(utilpath);
end

function [ok, user]=Login
prompt={'Please enter your username'};
name='Login';
numlines=1;
defaultanswer={'lab'};
user=inputdlg(prompt,name,numlines,defaultanswer);
if isempty(user) %user pressed cancel
    fprintf('\nUser pressed cancel, goodbye.')
    ok=0;
else
    ok=1;
end
end

function InitCamera(hObject, handles)
global pref;
try
    vid = videoinput('dalsa', 1, pref.ccf);
    src = getselectedsource(vid);
    %imaqmem(1e12);
    vid.timeout=60;
    vid.LoggingMode = 'memory';
    triggerconfig(vid,pref.trigger_type{1},pref.trigger_type{2},pref.trigger_type{3});
    vid.FramesPerTrigger = 1;
    Message('Camera initialized successfully', handles)
catch
    Message('Camera initialization failed', handles)
    questdlg('Failed to initialize camera! Running in simulation mode.',     'Camera Init Failure!', 'OK', 'Cancel', 'OK');
    vid=[];
end
handles.vid=vid;
guidata(hObject, handles)
end

function InitSoundOut(handles)
global pref;
userdata=get(handles.figure1, 'userdata');

% Start initializing while we get some information
InitializePsychSound(0);

% Get soundcard params & save
SoundFs         = pref.fs;
OutputDeviceID  = pref.dev_id;
numChan         = pref.n_chan;
buffSize        = pref.buff_size;
reqlatencyclass = 1;
handles.SoundFs.String = num2str(SoundFs);
set(handles.OutputDeviceID,'Value',OutputDeviceID+1);

%stop and close if already inited before
try
    PsychPortAudio('Stop', userdata.paOuthandle);
    PsychPortAudio('Close');
end

% Open psychportaudio, set run mode, and save handle.
try paOuthandle = PsychPortAudio('Open', OutputDeviceID, 1, reqlatencyclass, SoundFs, numChan, buffSize);
    PsychPortAudio('RunMode', paOuthandle, pref.runMode);
    
    userdata.paOuthandle=paOuthandle;
    set(handles.figure1, 'userdata', userdata);
    
    Message(sprintf('using Sound Output device %d', OutputDeviceID), handles)
    Message('Initialized Sound Output', handles) 
catch
    Message(sprintf('Error: could not open Output Device'), handles);
end
end     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Stimulus Handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function InitStim(stimuli, handles)
global pref
Message('', handles)
Message('Initializing new stimulus protocol', handles)

%load calibration into userdata
cd(pref.calibration)
try
    cal=load('calibration.mat');
    Message('Loaded speaker calibration data', handles)
    str=sprintf('Last calibration was flat to within std= +- %.1f dB (range %.1f - %.1f dB)', std(cal.DB), min(cal.DB), max(cal.DB));
    Message(str, handles)
catch
    Message('Error: could not load speaker calibration data. Tones will be uncalibrated', handles)
    cal=[];
end

% Compute information for frame triggers
Fs           = str2num(get(handles.SoundFs, 'string'));
triglength   = round(Fs/1000); %1 ms trigger
FPS          = pref.fps;
ttl_interval = 1000/FPS;
ttl_int_samp = round(ttl_interval*Fs/1000); %ttl_interval in samples

% Set stimulus info in GUI
str=stimuli(1).param.name;
Message(str, handles)
str=stimuli(1).param.description;
Message(str, handles)

% Make stimulus & triggers concurrently so they can be aligned.
toneseries=[];
triggerseries=[];
ledseries=[];
total_frames = 0;
for n=2:length(stimuli)
    %Make Tone
    typeidx = strcmp(pref.stimulitypes(:,1),stimuli(n).type);
    typefcn = pref.stimulitypes(typeidx,2);
    typefcn = typefcn{:};
    stimuli(n).param = CalibrateSound(stimuli(n).param, stimuli(n).type, cal, handles);
    sample = feval(typefcn,stimuli(n).param,Fs);
    
    % Make LED & Frame triggers
    stim_length  = length(sample);
    isi          = stimuli(n).param.next;
    silence      = zeros(1, round(Fs*.001*(isi)));
    total_length = stim_length + length(silence);
    stim_frames  = ceil(total_length/ttl_int_samp); 
    new_length   = stim_frames*ttl_int_samp;
    
    % Since we want too many rather than too few frames, pad end of silence
    if new_length>total_length
        diff = new_length-total_length;
        silence((end+1):(end+diff)) = zeros(diff,1);
        total_length = new_length;
    end

    sample_triggers = zeros(total_length,1);
    for i=0:(stim_frames-1)
        sample_triggers((i*ttl_int_samp)+1:((i*ttl_int_samp)+triglength)) = [0.25,.5*ones(size(2:triglength))];
    end

    % Append
    sample_led    = .5*ones(size(sample_triggers'));
    toneseries    = [toneseries, sample, silence];
    triggerseries = [triggerseries, sample_triggers'];
    ledseries     = [ledseries, sample_led];
    total_frames  = total_frames + stim_frames;
end

% For now, not using green channel so pad with zeros
greenchan = zeros(size(toneseries));

% Combine all vectors to tone matrix
tone=zeros(length(toneseries),4); %iti is implemented as silence after tone
tone(:,1)=toneseries;
tone(:,2)=triggerseries;
tone(:,3)=ledseries;
tone(:,4)=greenchan;

% Gather summary information on stimulus
timestamp             = datestr(now);
serieslength          = length(toneseries);
nreps                 = str2num(get(handles.numreps, 'string'));
duration              = serieslength/Fs;
total_duration        = nreps*serieslength/Fs;
series_periodicity    = Fs/serieslength;
total_duration_frames = nreps*total_frames;

% Print summary info to GUI
str=sprintf('Video frame interval: %d samples = %.4f ms = %.4f fps', ttl_int_samp, 1000*ttl_int_samp/Fs,Fs/ttl_int_samp);
Message(str, handles)
str=sprintf('This stimulus will collect %d frames', total_frames);
Message(str, handles)

%update StimOutputTable
data=[duration; series_periodicity; length(stimuli)-1; total_duration_frames; total_duration;total_duration/60 ];
set(findobj('tag', 'StimOutputTable'), 'Data', data)
stimparams.series_period_sec=duration;
stimparams.series_period_frames=total_frames;
stimparams.series_periodicity=series_periodicity;
stimparams.FPS=FPS;
%save stimparams
userdata=get(handles.figure1, 'userdata');
userdata.stimuli=stimuli;
userdata.stimparams=stimparams;
userdata.tone=tone;
set(handles.figure1, 'userdata', userdata);
end

function Go_Callback(hObject, eventdata, handles)
global pref
% hObject    handle to Go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and User data (see GUIDATA)

userdata=get(handles.figure1, 'userdata');

%reset abort flag
userdata.abort=0;
set(handles.figure1, 'userdata', userdata);

%load stimuli
if isfield(userdata, 'tone')
    % Do nothing because stimulus is already inited
elseif isfield(userdata, 'stimuli')
    % Stim file is somehow loaded without being inited
    InitStim(userdata.stimuli, handles);
else
    error('no stimuli selected yet')
end

vid=handles.vid;
if isrunning(vid)
    stop(vid);
end

if isempty(vid) Message('\nsimulation mode', handles);return; end

% Set Camera Options
vid.FramesPerTrigger = 1;
vid.TriggerRepeat = Inf;
triggerconfig(vid, pref.trigger_type{1},pref.trigger_type{2},pref.trigger_type{3});


%Notify user of start
data=get(findobj('tag', 'StimOutputTable'), 'Data');
filenum=userdata.filenum;
nframes=data(4);
totaldurationsecs=data(5);
%vid.timeout=totaldurationsecs+10;
Message(sprintf('\nInitializing camera to collect %d frames (%.1f s)', nframes, totaldurationsecs), handles)

% Prepare live previewing figure
% See http://www.mathworks.com/help/imaq/examples/video-display-with-live-histogram.html

prevfig = figure('Visible','off');
image_res = fliplr(vid.VideoResolution);
subplot(1,2,1);
histplot = imshow(zeros(image_res));
axis image;
setappdata(histplot,'UpdatePreviewWindowFcn',@update_livehistogram_display);
preview(vid,histplot);
handles.prevfig = prevfig;

% Make videowriter object & attach to the video input object
cd(userdata.datadir)
filenum=userdata.filenum;
fn=sprintf('M-%d', filenum);
fnev=sprintf('M-%d-events', filenum);
fnud=sprintf('M-%d-stimparams',filenum);
vwrite = VideoWriter(fn,'Archival');
vwrite.FrameRate = pref.fps;
vid.LoggingMode = 'disk';
vid.DiskLogger = vwrite;


%Start video object
start(vid);
preview(vid)
if isrunning(vid)
    Message('imaq started, waiting for triggers...', handles)
    Message('starting stimulus...', handles)
else
    Message('Video object was not started! Pausing 3s and trying again', handles)
    pause(3);
    if isrunning(vid)
        Message('imaq started, waiting for triggers...', handles)
    else
        Message('Object could not be started, not starting stimulus.', handles)
        return
    end
end

PlaySound(handles)

% Wait until we're done
status = PsychPortAudio('GetStatus',handles.figure1.UserData.paOuthandle);
while status.Active == 1
    status = PsychPortAudio('GetStatus',handles.figure1.UserData.paOuthandle);
    drawnow;
end

Message('Stimuli complete! Waiting for disk writer to catch up...', handles)
while (vid.FramesAcquired ~= vid.DiskLoggerFrameCount)
    pause(1)
    Message(sprintf('Total frames: %d, Logger Frames: %d',vid.FramesAcquired,vid.DiskLoggerFrameCount),handles)
    drawnow
end

Message('Video saved, writing event & stimulus information',handles)
drawnow
evts = vid.EventLog;
save(fnev, 'evts','-v7.3');
save(fnud, 'userdata','-v7.3');
userdata.filenum = filenum+1;
set(handles.figure1, 'userdata', userdata);

%Delete preview window & clean up
try
delete(handles.prevfig);
end

stop(vid);
%Save_Frames(handles);
Message('Success!',handles)
drawnow

end

function PlaySound(handles)
userdata=get(handles.figure1, 'userdata');
tone=userdata.tone;
numreps=str2num(get(handles.numreps, 'string'));
for n=1:numreps
    if ~userdata.abort
        paOuthandle=userdata.paOuthandle;
        PsychPortAudio('FillBuffer', paOuthandle, tone'); % fill buffer now, start in PlaySound
        nreps=1;
        when=0; %use this to start immediately
        waitForStart=0;
        PsychPortAudio('Start', paOuthandle,nreps,when,waitForStart);
        Message(sprintf('rep %d', n), handles)
    end 
end
handles.figure1.UserData.playing = 0;
end

function stimparams=CalibrateSound(stimparams, stimtype, cal, handles)
if ~isempty(cal) %it will be empty if Init failed to load calibration
    if strcmp(stimtype, '2tone') %special case since 2tone has both a frequency and a probefreq
        try
            findex=find(cal.logspacedfreqs<=stimparams.frequency, 1, 'last');
            atten=cal.atten(findex);
            stimparams.amplitude=stimparams.amplitude-atten;
           
            findex=find(cal.logspacedfreqs<=stimparams.probefreq, 1, 'last');
            atten=cal.atten(findex);
            stimparams.probeamp=stimparams.probeamp-atten;
           
            %Message('calibrated', handles)
        catch
            Message('NOT calibrated', handles)
        end
        
        
    elseif isfield(stimparams, 'frequency') %it has a freq and therefore is calibratable by frequency
        try
            findex=find(cal.logspacedfreqs<=stimparams.frequency, 1, 'last');
            atten=cal.atten(findex);
            switch stimtype
                case 'bintone'
                    Ratten=cal.Ratten(findex);
                    Latten=cal.Latten(findex);
                    stimparams.Ramplitude=stimparams.Ramplitude-Ratten;
                    stimparams.Lamplitude=stimparams.Lamplitude-Latten;
                otherwise
                    stimparams.amplitude=stimparams.amplitude-atten;
            end
            %Message('calibrated', handles)
        catch
            Message('NOT calibrated', handles)
        end
        
    else
        switch stimtype
            case {'clicktrain', 'whitenoise', 'amnoise'} %stimuli that consist of white noise
                try
                    findex=find(cal.logspacedfreqs==-1); %freq of -1 indicates white noise
                    atten=cal.atten(findex);
                    switch stimtype
                        case 'binwhitenoise'
                            Ratten=cal.Ratten(findex);
                            Latten=cal.Latten(findex);
                            stimparams.Ramplitude=stimparams.Ramplitude-Ratten;
                            stimparams.Lamplitude=stimparams.Lamplitude-Latten;
                        otherwise
                            stimparams.amplitude=stimparams.amplitude-atten;
                    end
                    %Message(sprintf('calibrated'), handles)
                catch
                    Message('NOT calibrated', handles);
                end
            case {'fmtone'} %stimuli that have a carrier frequency
                try
                    findex=find(cal.logspacedfreqs<=stimparams.carrier_frequency, 1, 'last');
                    atten=cal.atten(findex);
                    stimparams.amplitude=stimparams.amplitude-atten;
                    %Message('calibrated', handles)
                catch
                    Message('NOT calibrated', handles);
                end
            case {'noise'} %narrow-band noise stimuli (use center frequency calibration)
                try
                    findex=find(cal.logspacedfreqs<=stimparams.center_frequency, 1, 'last');
                    atten=cal.atten(findex);
                    stimparams.amplitude=stimparams.amplitude-atten;
                    %Message('calibrated', handles)
                catch
                    Message('NOT calibrated', handles)
                end
            case {'GPIAS'} %startle pulse (use whitenoise calibration)
                %plus narrow-band noise (use center frequency calibration)
                try
                    findex=find(cal.logspacedfreqs<=stimparams.center_frequency, 1, 'last');
                    atten=cal.atten(findex);
                    stimparams.amplitude=stimparams.amplitude-atten;
                    findex2=find(cal.logspacedfreqs==-1); %freq of -1 indicates white noise
                    atten=cal.atten(findex2);
                    stimparams.pulseamp=stimparams.pulseamp-atten;

                    %Message('calibrated', handles)
                catch
                    Message('NOT calibrated', handles)
                end
                case {'ASR'} %startle pulse (use whitenoise calibration)
                %plus narrow-band noise pulse (use center frequency calibration)
                
                try % there is a bug here!!! mak 26June2012
                    % findex finds the index before the index desired.
                    % E.g., for 4000 Hz it gives 3364 Hz, although it's
                    % only a difference of 4 dB
                    findex=find(cal.logspacedfreqs<=stimparams.prepulsefreq, 1, 'last'); 
                    atten=cal.atten(findex);
                    stimparams.prepulseamp=stimparams.prepulseamp-atten;
                    findex2=find(cal.logspacedfreqs==-1); %freq of -1 indicates white noise
                    atten=cal.atten(findex2);
                    stimparams.pulseamp=stimparams.pulseamp-atten;

                    %Message('calibrated', handles)
                catch
                    Message('NOT calibrated', handles)
                end
        end
    end
end
end

function AddStimtoList(stim, handles)
userdata=get(handles.figure1, 'userdata');
h=handles.CurrentStimulusProtocol;
str=get(h, 'String');
if isfield(userdata, 'StimList')
StimList=userdata.StimList;
else
    StimList=[];
end

if ~iscell(str) s{1}=str; str=s;end
numfiles=length(str);
current=stim.stimuli(1).param.name;
newstr{1}=current;
newStimList(1)=stim;
history_length=str2num(get(handles.StimHistoryLength, 'string'));
for n=1:(min(history_length-1, length(StimList)))
    if isempty(str{n}) | strcmp(str{n}, ' ') | strcmp(str{n}, '') %skip
    else
        newstr{n+1}=str{n};
        if isempty(StimList)
        else
            newStimList(end+1) = StimList(n);
        end
    end
end
set(h, 'string', newstr, 'value', 1)
pathname=fileparts(which(mfilename));
cd(pathname)
save RecentStimulusProtocols  newstr newStimList
userdata.StimList=newStimList;
set(handles.figure1, 'userdata', userdata);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Data & Utility Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CheckLightLevels_Callback(hObject, eventdata, handles)
% hObject    handle to CheckLightLevels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and User data (see GUIDATA)

% Make image window
BG=get(hObject, 'backgroundcolor');
Str=get(hObject, 'string');
set(hObject, 'backgroundcolor', 'r', 'string', 'one moment...')
close(findobj('Type', 'figure', 'Tag', 'Light Level Histogram'));
close(findobj('Type', 'figure', 'Tag', 'Image'));

% Take and average ten frames
vid=handles.vid;
if isempty(vid) fprintf('\nsimulation mode');return; end
vid.FramesPerTrigger=10;
triggerconfig(vid, 'immediate');
start(vid)
m = getdata(vid, 10);
stop(vid)
m=squeeze(m);
m=mean(m,3);

% Draw window
figure;
set(gcf, 'pos', [44   271   560   420])
imagesc(m)
colormap(gray)
title('Image')
set(gcf, 'Tag', 'Image');
fig=figure;
set(gcf, 'pos', [44   771   560   250])
hist(reshape(m, 1, prod(size(m))), 1000);
title('Light Level Histogram')
set(gcf, 'Tag', 'Light Level Histogram');
set(hObject, 'backgroundcolor', BG, 'string', Str)
end

function CreateDataDir(handles)
global pref;
warning off MATLAB:MKDIR:DirectoryExists
dataroot=pref.data;
cd(dataroot)

expdate=datestr(now, 'mmddyy');
if ~exist(expdate, 'dir')
    mkdir(expdate)
end
cd(expdate)


sess_idx=1;
session=sprintf('%03d',sess_idx);
while exist(session, 'dir')
    sess_idx=sess_idx+1;
    session=sprintf('%03d',sess_idx);
end
mkdir(session)
cd(session)
h=findobj('Tag', 'DataPath');
set(h, 'string', pwd)

%store current data dir in userdata
userdata=get(handles.figure1, 'userdata');
userdata.datadir=pwd;
userdata.filenum=0;
set(handles.figure1, 'userdata', userdata);
end

function Save_Frames(handles)
userdata=get(handles.figure1, 'userdata');
vid=handles.vid;

%Stop waiting for triggers
if isrunning(vid)
    stop(vid);
else
    Message('Video object already stopped when we went to save',handles);
end

%Get frames from memory buffer
nframes = vid.FramesAvailable;
if nframes == 0
    Message('No frames were recorded, not saving images.', handles)
else
    try
        nframes = vid.FramesAvailable;
        filenum=userdata.filenum;
        fn=sprintf('M-%d.mat', filenum);
        fnts=sprintf('M-%d-timestamps.mat', filenum);
        Message(sprintf('Retrieved %d Frames. Saving raw video data as %s',nframes,fn), handles)
        drawnow;
        [M,timestamps] = getdata(vid, nframes);
        cd(userdata.datadir)
        save(fn, 'M','-7.3');
        save(fnts, 'timestamps');
        Message('done', handles)
        userdata.filenum = filenum+1;
    catch
        Message('Error getting frames!',handles);
        %M = getdata(vid, get(vid, 'FramesAvailable'));
    end
end

%check the AnalyzeWhenDone checkbox
% Hint: get(hObject,'Value') returns toggle state of AnalyzeWhenDone
hAnalyzeWhenDone=findobj('Tag', 'AnalyzeWhenDone');
if get(hAnalyzeWhenDone,'Value')
    %     imanal2(pwd)
    fprintf('\nAnalyzing...')
    fft_mem(pwd, M);
end
%write video data to file for possible later analysis
%cd to current data dir
cd(userdata.datadir)
set(handles.figure1, 'userdata', userdata);

% Write stimulus info to file
save stimparams userdata


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Secondary Callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Abort_Callback(hObject, eventdata, handles)
% hObject    handle to Abort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(handles.figure1, 'userdata');
%get(hObject,'value')
try
    Message('aborting run...', handles)
    stop(handles.vid);
    set(handles.vid, 'timeout', 1);
    PsychPortAudio('Stop', userdata.paOuthandle);
    userdata.abort=1;
    set(handles.figure1, 'userdata', userdata);
%     set(hObject, 'backgroundcolor', 'r')
catch
    Message('failure during abort', handles)

end
end

function AnalyzeCurrentDir_Callback(hObject, eventdata, handles)
% hObject    handle to AnalyzeCurrentDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%imanal2(pwd)
fft_mem(pwd);
end

function AnalyzeWhenDone_Callback(hObject, eventdata, handles)
% hObject    handle to AnalyzeWhenDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AnalyzeWhenDone
end

function ChangeDir_Callback(hObject, eventdata, handles)
% hObject    handle to ChangeDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newpath = uigetdir('..', 'Choose new directory in which to save data')
if newpath
    cd(newpath)
    h=findobj('Tag', 'DataPath');
    set(h, 'string', pwd)
    
    %store current data dir in userdata
    userdata=get(handles.figure1, 'userdata');
    userdata.datadir=pwd;
    set(handles.figure1, 'userdata', userdata);

end
end

function ClearMessages_Callback(hObject, eventdata, handles)
% hObject    handle to ClearMessages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=handles.Message;
set(h, 'String', '');
end

function ClearStimHistory_Callback(hObject, eventdata, handles)
% hObject    handle to ClearStimHistory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userdata=get(handles.figure1, 'userdata');
h=handles.CurrentStimulusProtocol;
set(h, 'String', []);
pathname=fileparts(which(mfilename));
cd(pathname)
delete('RecentStimulusProtocols.mat')  
userdata.StimList=[];
set(handles.figure1, 'userdata', userdata);


%myabe this might help close all devices if getting unavailable errors
%for n=1:30;try(PsychPortAudio('Stop', n)), end, end
%for n=1:30;try(PsychPortAudio('Close', n)), end, end
end

function CurrentStimulusProtocol_Callback(hObject, eventdata, handles)
% hObject    handle to label14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns label14 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from label14
%set(hObject, 'string', '', 'value', 1)
%set(hObject, 'value', 1)
    userdata=get(handles.figure1, 'userdata');
    StimList=userdata.StimList;
    value=get(hObject, 'value');
%     pathname=StimList(value).pathname;
%     filename=StimList(value).filename;
stimuli=StimList(value).stimuli;
    
%     global pref
%     cd(pref.experhome)
%     cd('protocols')
%     cd('Tuning Curve protocols')
%     files=get(hObject, 'string');
%     value=get(hObject, 'value');
%     filename=files{value};
%     load(filename);
    InitStim(stimuli, handles);
end

function LoadStim_Callback(hObject, eventdata, handles)
% hObject    handle to LoadStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pref
cd(pref.protocols)
[filename, pathname]=uigetfile('*.mat', 'Choose Stimulus Protocol');
cd(pathname)
load(filename);
stim.pathname=pathname;
stim.filename=filename;
stim.stimuli=stimuli;
InitStim(stimuli, handles);
AddStimtoList(stim, handles)
end

function NewDir_Callback(hObject, eventdata, handles)
% hObject    handle to NewDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CreateDataDir(handles)
end

function numreps_Callback(hObject, eventdata, handles)
% hObject    handle to numreps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numreps as text
%        str2double(get(hObject,'String')) returns contents of numreps as a double
userdata=get(handles.figure1, 'userdata');
if isfield(userdata, 'stimuli')
InitStim(userdata.stimuli, handles);
end
end

function OutputDeviceID_Callback(hObject, eventdata, handles)
% hObject    handle to OutputDeviceID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns OutputDeviceID contents as cell array
%        contents{get(hObject,'Value')} returns selected item from OutputDeviceID
OutputDeviceID = get(handles.OutputDeviceID, 'Value')-1;
userdata=get(handles.figure1, 'userdata');
userdata.OutputDeviceID=OutputDeviceID;
set(handles.figure1, 'userdata', userdata);
InitSoundOut(handles)
end

function Reset_Callback(hObject, eventdata, handles)
% hObject    handle to Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(handles.figure1, 'userdata');
try
    vid=handles.vid;
stop(vid)
clear vid
InitCamera(hObject, handles);
InitSoundOut(handles)
Message('Reset successful', handles)
userdata.abort=0;
set(handles.figure1, 'userdata', userdata);

catch
    Message('failure during reset', handles)

end
end

function SetRoot_Callback(hObject, eventdata, handles)
% hObject    handle to SetRoot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and User data (see GUIDATA)
end

function SoundFs_Callback(hObject, eventdata, handles)
% hObject    handle to SoundFs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SoundFs as text
%        str2double(get(hObject,'String')) returns contents of SoundFs as a double
end

function StimHistoryLength_Callback(hObject, eventdata, handles)
% hObject    handle to StimHistoryLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StimHistoryLength as text
%        str2double(get(hObject,'String')) returns contents of StimHistoryLength as a double
end

function User_Callback(hObject, eventdata, handles)
% hObject    handle to User (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and User data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of User as text
%        str2double(get(hObject,'String')) returns contents of User as a double
end

function ViewCurrentDir_Callback(hObject, eventdata, handles)
% hObject    handle to ViewCurrentDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imview(pwd)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Figure Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Message(string, handles)
h=handles.Message;
old_string=get(h, 'string');
if iscell(old_string)
    %old_string=old_string{:};end
    n= length(old_string);
    old_string{n+1}=string;
    set(h, 'String', old_string);
else
    new_string={old_string, string};
    set(h, 'String', new_string);
end
try
    jhEdit = findjobj(h);
    jEdit = jhEdit.getComponent(0).getComponent(0);
    jEdit.setCaretPosition(jEdit.getDocument.getLength);
end
end

function varargout = imageGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and User data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'vid')
    vid=handles.vid;
    delete(vid)
    clear vid
end
userdata=get(handles.figure1, 'userdata');

try
    % % Stop playback:
    PsychPortAudio('Stop', userdata.paOuthandle);
    % % Close the audio device:
    PsychPortAudio('Close'); % Excluding the handle shuts down the entire driver.
end
% Hint: delete(hObject) closes the figure
delete(hObject);
end

function CurrentStimulusProtocol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to label14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
pathname=fileparts(which(mfilename));
cd(pathname)
try
    load('RecentStimulusProtocols.mat')
    set(hObject, 'String', newstr, 'value', 1);
   % stimuli=newStimList(1).stimuli;
   % InitStim(stimuli, handles);
end
end

function numreps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numreps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function OutputDeviceID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OutputDeviceID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
InitializePsychSound(0);
% note: in the list of devices, the first one is device0
% (devs(1).DeviceIndex=0)
devs = PsychPortAudio('GetDevices');
for i = 1:length(devs)
    deviceString{i}=sprintf('%d: %s: %s', devs(i).DeviceIndex, devs(i).HostAudioAPIName, devs(i).DeviceName);
end
set(hObject, 'String', deviceString);
set(hObject, 'Value',29); %default Output Device ID, note devs is 0-indexed
end

function SoundFs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SoundFs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function StimHistoryLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StimHistoryLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function User_CreateFcn(hObject, eventdata, handles)
% hObject    handle to User (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


