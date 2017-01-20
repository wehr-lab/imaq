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

% Last Modified by GUIDE v2.5 13-Jan-2017 15:23:36

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
fs       = pref.fs;

%%% Compute derived values
% Generate list of frequencies
% 12 freqs per octave for now, we'll see if it matters
noctaves=log2(maxfreq/minfreq);
freq=minfreq*2.^([0:(1/20):noctaves]); 

% spectrogram values
nseg = dur/4; %number of segments w/ a 12.5ms window
segsamples = round(((dur/1000)*fs)/nseg); % How many samples in a segment
noverlap = round(segsamples); % 1/3 overlap
winsize = segsamples+noverlap; % window size. Don't ask me why it's computed this way
nfft = 10000;

% make bandpass filter
lowrolloff  = minfreq/(2^(1/3)); % hardcoding a major third for now (~.8, 1.25 freq ratio)
highrolloff = min(maxfreq*(2^(1/3)),fs/2); % the calculated rolloff or nyquist.
ast = 60; % Amount of attenuation in the stop bands (dB)
ap  = 1;  % Amount of ripple allowed in the passband
bpf_specs = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',lowrolloff,minfreq,maxfreq,highrolloff,ast,ap,ast,fs);
bpf = design(bpf_specs, 'butter');

for i = 1:niter
    %%% Synthesize sound
    handles.status_text.String = 'Generating Sound';
    drawnow
    snd = genripple(rand*am_vel,rand*fm_dens,rand*depth,dur,fs,freq);
    
    % Bandpass filter sound
    

    % Generate & plot synthesized spectrogram
    handles.status_text.String = 'Generating Spectrogram';
    axes(handles.synthesized);
    spectrogram(snd,winsize,noverlap,nfft,fs,'yaxis');
    %set(gca,'YScale','log');
    ylim([minfreq/1000,maxfreq/1000]);
    colorbar('off');

    % Generate & plot synthesized spectrum
    handles.status_text.String = 'Generating Spectrum';
    axes(handles.filter);
    %periodogram(snd,[],nfft,fs);
    %xlim([minfreq/1000,maxfreq/1000]);
    [Pxx,f] = pwelch(snd,winsize,noverlap,nfft,fs);
    plot(f,Pxx)
    xlim([minfreq,maxfreq])
    set(gca,'YScale','log')

    %%% Apply filter
    if i>1
        snd = snd';
        snd_out = zeros(length(snd),1);
        snd = [zeros(192-1,1); snd ];
        for n = 1:length(snd_out)
            % Convolve filter w/ input signal
            snd_out(n) = weights' * snd(n:n+192-1);
        end
        
        snd = snd_out';
        snd = ((snd)./(max(max(snd),abs(min(snd))))).*0.5;
    end
        
    
    %%% Play & Record sound
    handles.status_text.String = 'Playing Audio';
    drawnow
    pahandle=handles.audio;
    PsychPortAudio('FillBuffer', pahandle, snd); % fill playback buffer
    PsychPortAudio('GetAudioData', pahandle, dur/1000); % Preallocate a recording buffer
    
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
    % Normalize recorded sound like synthesized sound
    recorded = ((recorded)./(max(max(recorded),abs(min(recorded))))).*0.5;
    
    %%% Plot Measured sound
    handles.status_text.String = 'Generating Spectrogram';
    axes(handles.measured);
    spectrogram(recorded,winsize,noverlap,nfft,fs,'yaxis');
    %set(gca,'YScale','log');
    ylim([minfreq/1000,maxfreq/1000]);
    colorbar('off');
    
    handles.status_text.String = 'Generating Spectrum';
    axes(handles.filter);
    hold on
    %periodogram(snd,[],nfft,fs);
    %xlim([minfreq/1000,maxfreq/1000]);
    [Pxx,f] = pwelch(recorded,winsize,noverlap,nfft,fs);
    plot(f,Pxx,'g')
    
    %%% Compute filter
    if i == 1
        [signal_out, err, weights] = lms_filter(recorded, snd,[]);
    else
        [signal_out, err, weights] = lms_filter(recorded, snd, weights);
    end
    signal_out = ((signal_out)./(max(max(signal_out),abs(min(signal_out))))).*0.5;
    [Pxx,f] = pwelch(signal_out,winsize,noverlap,nfft,fs);
    plot(f,Pxx,'r')
    hold off
end

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

function [signal_out, err, weights] = lms_filter(signal_in, desired, weights)
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
FilterLength = 192;              % lowest we need to filter will be ~1kHz (192k/192)
mu = 125e-6;                     % Adaptation step size:

if isempty(weights)
    weights = zeros(FilterLength,1); % Filter coefficients:
end
        
% Pre-allocate output and error signals:
signal_out = zeros(nsamples,1);
err = zeros(nsamples,1);
       
%%% Compute filter
% Zero Pad Input Signal:
signal_in = [zeros(FilterLength-1,1); signal_in ];
for n = 1:nsamples
    
    % Convolve filter w/ input signal
    signal_out(n) = weights' * signal_in(n:n+FilterLength-1);
    
    % Update the filter coefficients:
    err(n) = desired(n) - signal_out(n) ;
    weights = weights + mu*err(n)*signal_in(n:n+FilterLength-1);
    
end

    

    


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



function freqsperoct_val_Callback(hObject, eventdata, handles)
% hObject    handle to freqsperoct_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freqsperoct_val as text
%        str2double(get(hObject,'String')) returns contents of freqsperoct_val as a double


% --- Executes during object creation, after setting all properties.
function freqsperoct_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freqsperoct_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end