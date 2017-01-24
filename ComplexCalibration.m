function varargout = ComplexCalibration(varargin)
% COMPLEXCALIBRATION MATLAB code for ComplexCalibration.fig
%      COMPLEXCALIBRATION, by itself, creates a new COMPLEXCALIBRATION or raises the existing
%      singleton*.
%
%      H = COMPLEXCALIBRATION returns the handle to a new COMPLEXCALIBRATION or the handle to
%      the existing singleton*.
%
%      COMPLEXCALIBRATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COMPLEXCALIBRATION.M with the given input arguments.
%
%      COMPLEXCALIBRATION('Property','Value',...) creates a new COMPLEXCALIBRATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ComplexCalibration_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ComplexCalibration_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ComplexCalibration

% Last Modified by GUIDE v2.5 20-Jan-2017 15:26:23

% Begin initialization code - DO NOT EDIT

%%%%%%%%%%%%%%%
%%%%
%%%%
%%%% ADD PREFS CALL
%%%% CHANGE UI FOR DEPTH TO 0-1
%%%% GET FS from prefs in GO CB
%%%%
%%%%%%%%%%%%%%%



gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ComplexCalibration_OpeningFcn, ...
                   'gui_OutputFcn',  @ComplexCalibration_OutputFcn, ...
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


% --- Executes just before ComplexCalibration is made visible.
function ComplexCalibration_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ComplexCalibration (see VARARGIN)

global pref
Prefs('nologin');

% Choose default command line output for ComplexCalibration
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Get & Set values for plots
minfreq = str2num(handles.minfreq_val.String);
maxfreq = str2num(handles.maxfreq_val.String);
maxdur  = str2num(handles.dur_val.String);

lims = [minfreq,maxfreq];
dur  = [0,maxdur];

handles.synthesized.YLim = lims;
handles.measured.YLim    = lims;
handles.filter.XLim      = lims;
handles.synthesized.XLim = dur;
handles.measured.XLim    = dur;

% Field for stored kernel
handles.kernel_weights = [];

% Running state indicator

InitSound(hObject,handles);


% --- Outputs from this function are returned to the command line.
function varargout = ComplexCalibration_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in go_button.
function go_button_Callback(hObject, eventdata, handles)
% hObject    handle to go_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%
%%% Overview
% 1: Collect necessary values
%%% Iterate...
% 2: Generate & play ripple sound
% 3: Compute inverse filter & plot
% 4: Apply inverse filter and repeat
%%%%%%%%%%%%

% Change button appearance
handles.go_button.String = 'STOP';
handles.go_button.BackgroundColor = [1,0,0];


%%%%% 1: Collect values
global pref
minfreq  = str2num(handles.minfreq_val.String);
maxfreq  = str2num(handles.maxfreq_val.String);
dur      = str2num(handles.dur_val.String);
am_vel   = str2num(handles.am_val.String);
fm_dens  = str2num(handles.fm_val.String);
depth    = str2num(handles.mod_val.String);
niter    = str2num(handles.niter_val.String);
converge = str2num(handles.converge_val.String);
k_size   = str2num(handles.kernel_val.String);
weights  = handles.kernel_weights;
fs       = pref.fs;

%%% Compute derived values


% spectrogram values
nseg = dur/4; %number of segments w/ a 12.5ms window
segsamples = round(((dur/1000)*fs)/nseg); % How many samples in a segment
noverlap = round(segsamples); % 1/3 overlap
winsize = segsamples+noverlap; % window size. Don't ask me why it's computed this way
nfft = 10000;

% make bandpass filter
% lowrolloff  = minfreq/(2^(1/3)); % hardcoding a major third for now (~.8, 1.25 freq ratio)
% highrolloff = min(maxfreq*(2^(1/3)),fs/2); % the calculated rolloff or nyquist.
% ast = 60; % Amount of attenuation in the stop bands (dB)
% ap  = 1;  % Amount of ripple allowed in the passband
% bpf_specs = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',lowrolloff,minfreq,maxfreq,highrolloff,ast,ap,ast,fs);
% bpf = design(bpf_specs, 'butter');
% handles.status_text.String = 'Making Bandpass Filter';
% bpf = designfilt('bandpassiir',...
%     'FilterOrder',128,...
%     'SampleRate',fs,...
%     'HalfPowerFrequency1',minfreq-100,...
%     'HalfPowerFrequency2',maxfreq+100,...
%     'DesignMethod','butter');

for i = 1:niter
    if handles.go_button.Value == 0
       handles.weights = weights;
       handles.go_button.String = 'GO';
       handles.go_button.BackgroundColor = [0.114,0.678,0.282];
       break
    end
    
    % Generate list of frequencies
    % Generate list and randomly jitter spacing so carriers get spread.
    jitterfactor = rand*0.1;
    noctaves=log2((maxfreq-maxfreq*jitterfactor)/(minfreq+minfreq*jitterfactor));
    freq=(minfreq+minfreq*jitterfactor)*2.^([0:(1/24):noctaves]); 
    
    handles.thisiter_val.String = num2str(i);
    
    %%% Synthesize sound
    handles.status_text.String = 'Generating Sound';
    drawnow
    snd = genripple((rand*am_vel)+1,(rand*fm_dens)+1,rand*depth,dur,fs,freq);
    
    % Bandpass filter sound
    %snd = filtfilt(bpf,snd); % zero-phase filtering

    % Generate & plot synthesized spectrogram
    %handles.status_text.String = 'Generating Spectrogram';
    axes(handles.synthesized);
    spectrogram(snd,winsize,noverlap,nfft,fs,'yaxis');
    %set(gca,'YScale','log');
    ylim([minfreq/1000,maxfreq/1000]);
    colorbar('off');

    % Generate & plot synthesized spectrum
%     handles.status_text.String = 'Generating Spectrum';
%     axes(handles.filter);
%     %periodogram(snd,[],nfft,fs);
%     %xlim([minfreq/1000,maxfreq/1000]);
%     [Pxx,f] = pwelch(snd,winsize,noverlap,nfft,fs);
%     plot(f,Pxx)
%     xlim([minfreq,maxfreq])
%     set(gca,'YScale','log')

    %%% Apply filter
%     if i>1
%         snd = snd';
%         snd_out = zeros(length(snd),1);
%         snd = [zeros(k_size-1,1); snd ];
%         for n = 1:length(snd_out)
%             % Convolve filter w/ input signal
%             snd_out(n) = weights' * snd(n:n+192-1);
%         end
%         
%         snd = snd_out';
%         snd = ((snd)./(max(max(snd),abs(min(snd))))).*0.5;
%     end

     % weights are b coeffs in FIR. do zero-phase filtering
     %if ~isempty(weights)
     if i>100
         snd_out = filtfilt(weights,1,snd);
         % Normalize filtered sound
         powers = sqrt(bandpower([snd',snd_out'],fs,[minfreq,maxfreq]));
         snd_power = powers(1);
         out_power = powers(2);
         snd_out = snd_out*(snd_power/out_power);
     else
         snd_out = snd;
     end
     
    
    %%% Play & Record sound
    handles.status_text.String = 'Playing Audio';
    drawnow
    recorded = PlaySound(handles,snd_out,dur);
    
    %%% Plot Measured sound
    %handles.status_text.String = 'Generating Spectrogram';
    axes(handles.measured);
    spectrogram(recorded,winsize,noverlap,nfft,fs,'yaxis');
    %set(gca,'YScale','log');
    ylim([minfreq/1000,maxfreq/1000]);
    colorbar('off');
    
    %Normalize recorded sound
    powers = sqrt(bandpower([snd',recorded'],fs,[minfreq,maxfreq]));
    out_power = powers(1);
    rec_power = powers(2);
    snd = snd*(rec_power/out_power);
    
    % Get amplitude scaling
    % TBD
    
     %%% Compute filter
    handles.status_text.String = 'Computing Filter';
    drawnow
    if i == 1
        [signal_out, err, weights] = lms_filter(recorded, snd,[],k_size);
    else
        [signal_out, err, weights] = lms_filter(recorded, snd, weights,k_size);
    end   

%     % a lil test
%     if i == 1
%         [signal_out, err, weights] = lms_filter(snd, recorded,[],k_size);
%     else
%         [signal_out, err, weights] = lms_filter(snd, recorded, weights,k_size);
%     end   
    
    powers = sqrt(bandpower([signal_out,recorded'],fs,[minfreq,maxfreq]));
    out_power = powers(1);
    rec_power = powers(2);
    signal_out = signal_out*(rec_power/out_power);
    

    
    handles.status_text.String = 'Generating Spectrum';
    axes(handles.filter);
    cla
    hold on
    %periodogram(snd,[],nfft,fs);
    xlim([minfreq-100,maxfreq+100]);
    [Pxx_snd,f] = pwelch(snd,winsize,noverlap,nfft,fs);
    plot(f,Pxx_snd,'b')
    [Pxx_sndout,f] = pwelch(snd_out,winsize,noverlap,nfft,fs);
    plot(f,Pxx_sndout,'b--')
    [Pxx_rec,f] = pwelch(recorded,winsize,noverlap,nfft,fs);
    plot(f,Pxx_rec,'r')
    xlim([minfreq,maxfreq])
    set(gca,'YScale','log')
    [Pxx_filt,f] = pwelch(signal_out,winsize,noverlap,nfft,fs);
    plot(f,Pxx_filt,'g')
    hold off
end
% Save weights
weight_string = [pref.calibration,'calibration_kernel_',...
    datestr(now,'mmddyy-HHMM'),...
    '_F',num2str(minfreq),'-',num2str(maxfreq),...
    '_S',num2str(k_size),...
    '_N',num2str(niter),...
    '.mat'];
cal.weights = weights;
cal.minfreq = minfreq;
cal.maxfreq = maxfreq;
cal.power   = sqrt(bandpower(snd_out',fs,[minfreq,maxfreq]));
save(weight_string,'cal');
handles.kernel_weights = weights;
handles.status_text.String = sprintf('Finished %d Iterations',niter);
handles.go_button.String = 'GO';
handles.go_button.BackgroundColor = [0.114,0.678,0.282];
handles.go_button.Value = 0;



% Make sure to save weights!!!
    
    

   
   
%end
%%%% remember to send list of freqs generated from vals to gen_ripple







% --- Executes on button press in save_but.
function save_but_Callback(hObject, eventdata, handles)
% hObject    handle to save_but (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% Functions used during "Go"
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function snd = genripple(vel,dens,md,durrip,Fs,freq)
% Generate ripple sound, modified from "pa_genripple3" 
% in http://www.neural-code.com/index.php/panda

nFreq	= length(freq);
FreqNr	= 0:1:nFreq-1;
F0		= freq(1);
Oct		= log2(freq./F0);          % octaves above F0
Phi		= pi - 2*pi*rand(1,nFreq); % random phase
Phi(1)	= pi/2;                    % set first to 0.5*pi

nTime   = round((durrip/1000)*Fs); % Samples for Rippled Noise
time	= ((1:nTime)-1)/Fs;        % Time (sec)

% Generating the ripple
T		= 2*pi*vel*time;
F		= 2*pi*dens*Oct;
[T,F]	= meshgrid(T,F);
A		= 1+md*sin(T+F);
A		= A';

% Modulate carrier
snd		= 0; % 'initialization'
for ii = 1:nFreq
	rip		= A(:,ii)'.*sin(2*pi*freq(ii) .* time + Phi(ii));
	snd		= snd+rip;
end

% Normalization so max amp is +/- 0.5 w/o distorting sound
snd = ((snd)./(max(max(snd),abs(min(snd))))).*0.5;

function [signal_out, err, weights] = lms_filter(signal_in, desired, weights,k_size)
% LMS filter adapted from lms_01 in: 
% https://www.mathworks.com/help/simulink/ug/tutorial-integrating-matlab-code-with-a-simulink-model-for-filtering-an-audio-signal.html

nsamples = length(signal_in);

% If we get a row vector, make column vector
if size(signal_in,1)<size(signal_in,2)
    signal_in = signal_in';
end

if size(desired,1)<size(desired,2)
    desired = desired';
end

% Make sure the two input signals have the same length:
if any(size(desired) ~= size(signal_in))
    error([ 'The length of input argument ''desired'' is ' ...
        ' different from the length of ''signal_in''.' ]);
end
        
%%% Parameters
FilterLength = k_size;              % lowest we need to filter will be ~1kHz (192k/192)
mu = double(100e-7);                     % Adaptation step size:

if isempty(weights)
    weights = double(zeros(FilterLength,1)); % Filter coefficients:
else
    weights = double(weights);
end
        
% Pre-allocate output and error signals:
signal_out = zeros(nsamples,1);
err = double(zeros(nsamples,1));
desired = double(desired);
       
%%% Compute filter
% Zero Pad Input Signal:
signal_in = double([zeros(FilterLength-1,1); signal_in ]);
for n = 1:nsamples
    
    % Convolve filter w/ input signal
    signal_out(n) = weights' * signal_in(n:n+FilterLength-1);
    
    % Update the filter coefficients:
    err(n) =  desired(n) - signal_out(n) ;
%     if err(n)>0 && signal_out(n)>0
%         weights = weights - mu*err(n)*signal_in(n:n+FilterLength-1);

    weights = weights + mu*err(n)*signal_in(n:n+FilterLength-1); % Use VPA to not have rounding errors

%     elseif err(n)<=0 && signal_out(n)>0
%         weights = weights + mu*err(n)*signal_in(n:n+FilterLength-1);
%     elseif err(n)<=0 && signal_out(n)<=0
%         weights = weights - mu*err(n)*signal_in(n:n+FilterLength-1);
%     end
    
end

function recorded = PlaySound(handles,sound,dur)
    pahandle=handles.audio;
    
    % Have to pad with silence due to latency of soundcard, we
    % recover the real recording later.
    sound_pad = [zeros(48000,1)',sound,zeros(48000,1)'];
    
    PsychPortAudio('FillBuffer', pahandle, sound_pad); % fill playback buffer
    PsychPortAudio('GetAudioData', pahandle, (dur*1.5)/1000); % Preallocate a recording buffer
    
    nreps=1;when=0;waitForStart=0;
    PsychPortAudio('Start', pahandle,nreps,when,waitForStart);

    % Stop capture:
    waitForEndOfPlayback=1; %1: Wait until playback finished to stop
    PsychPortAudio('Stop', pahandle, waitForEndOfPlayback);

    % Retrieve audio data
    [recorded, ~, overflow, ~] = PsychPortAudio('GetAudioData', pahandle);
    if overflow>0
       handles.status_text.String = 'Warning! Overflow in recording!';
       drawnow
    end  
    
    % Find actual start of recording
    % mic is rl quiet before actual start, no sample exceeds 0.01, but first
    % sample of recording usually does.
    % We use diff because the mic can sometimes float above/below zero, a
    % strong departure is a better indicator than an absolute value
    startind = find(abs(diff(recorded))>0.008,1)+1; % Want one sample after so filter is causal
    recorded = recorded(startind:startind+length(sound)-1);

    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% Other accessory functions
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function InitSound(hObject,handles)
global pref
% InitializePsychSound(0); In imaq this is done by Prefs
%InitializeInputCh(handles) This gets done in InitParams
%InitializeOutputCh(handles)
SoundFs         = pref.fs;
DeviceID        = pref.dev_id;
buffSize        = pref.buff_size;
reqlatencyclass = 1;      % Try for lowest latency
numChan         = [1,1];  % Override prefs, One output, 1 input
runMode         = 1;      % Override prefs, leave soundcard on
mode            = 3;      % Full duplex: simultaneous capture and playback

% Check if Psychsound has been initialized (should be in Prefs)
try
    PsychPortAudio('GetDevices');
catch
    InitializePsychSound(1);
end

try pahandle = PsychPortAudio('Open', DeviceID, mode, reqlatencyclass, SoundFs, numChan, buffSize);
    PsychPortAudio('RunMode', pahandle, runMode);
    % Save handle of audio object
    handles.audio = pahandle;
    handles.status_text.String = 'Initialized sound successfully';
    guidata(hObject,handles)
catch
    handles.status_text.String = 'Could not load audio device!';
    handles.go_button.Enable = 'off';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% Useless shit that no one cares about
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function am_val_Callback(hObject, eventdata, handles)

function am_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to am_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dur_val_Callback(hObject, eventdata, handles)
% hObject    handle to dur_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

maxdur = str2num(get(hObject,'String'));
handles.synthesized.XLim(2) = maxdur;
handles.measured.XLim(2)    = maxdur;
pahandle = handles.audio;
PsychPortAudio('GetAudioData', pahandle, maxdur/1000);

function dur_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dur_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fm_val_Callback(hObject, eventdata, handles)

function fm_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fm_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxfreq_val_Callback(hObject, eventdata, handles)
% hObject    handle to maxfreq_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

maxfreq = str2num(get(hObject,'String'));
handles.synthesized.YLim(2) = maxfreq;
handles.measured.YLim(2)    = maxfreq;
handles.filter.XLim(2)      = maxfreq;

function maxfreq_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxfreq_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minfreq_val_Callback(hObject, eventdata, handles)
% hObject    handle to minfreq_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

minfreq = str2num(get(hObject,'String'));
handles.synthesized.YLim(1) = minfreq;
handles.measured.YLim(1)    = minfreq;
handles.filter.XLim(1)      = minfreq;

function minfreq_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minfreq_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mod_val_Callback(hObject, eventdata, handles)

function mod_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mod_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function niter_val_Callback(hObject, eventdata, handles)

function niter_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to niter_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function converge_val_Callback(hObject, eventdata, handles)

function converge_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to converge_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uitoolbar1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitoolbar1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
PsychPortAudio('Stop', handles.audio);
PsychPortAudio('Close', handles.audio);
end
% Hint: delete(hObject) closes the figure
delete(hObject);



function kernel_val_Callback(hObject, eventdata, handles)
% hObject    handle to kernel_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kernel_val as text
%        str2double(get(hObject,'String')) returns contents of kernel_val as a double


% --- Executes during object creation, after setting all properties.
function kernel_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kernel_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in reset_but.
function reset_but_Callback(hObject, eventdata, handles)
% hObject    handle to reset_but (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.kernel_weights = [];


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over go_button.
%function go_button_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to go_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in go_button.

% hObject    handle to go_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of go_button
