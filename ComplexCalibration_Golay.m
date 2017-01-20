function varargout = ComplexCalibration_Golay(varargin)
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

% Last Modified by GUIDE v2.5 13-Jan-2017 17:56:37

% Begin initialization code - DO NOT EDIT

%%%%%%%%%%%%%%%
%%%%
%%%% TODO:
%%%% 
%%%% 
%%%% 
%%%%
%%%%%%%%%%%%%%%



gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ComplexCalibration_Golay_OpeningFcn, ...
                   'gui_OutputFcn',  @ComplexCalibration_Golay_OutputFcn, ...
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
function ComplexCalibration_Golay_OpeningFcn(hObject, eventdata, handles, varargin)
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
minfreq = str2num(handles.lowfreq_val.String);
maxfreq = str2num(handles.highfreq_val.String);
lims = [minfreq,maxfreq];
handles.frequency.XLim      = lims;

InitSound(hObject,handles);


% --- Outputs from this function are returned to the command line.
function varargout = ComplexCalibration_Golay_OutputFcn(hObject, eventdata, handles) 
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
% 3: Compute inverse frequency & plot
% 4: Apply inverse frequency and repeat
%%%%%%%%%%%%

%%%%% 1: Collect values
global pref
order       = str2num(handles.order_val.String);
lowrolloff  = str2num(handles.lowrolloff_val.String);
lowfreq     = str2num(handles.lowfreq_val.String);
highfreq    = str2num(handles.highfreq_val.String);
highrolloff = str2num(handles.highrolloff_val.String);
rolloffamp  = str2num(handles.rolloffamp_val.String);
niter       = str2num(handles.niter_val.String);
kernelsize  = str2num(handles.kernel_val.String);
fs          = pref.fs;

% Calculate amplitude from previous calibration, if loaded
try
    noise_atten = pref.cal.atten(1);
    amplitude = pref.maxSPL-noise_atten;
    amplitude = 1*(10.^((amplitude-pref.maxSPL)/20));
catch
    amplitude = 0.5;
end

%for i = 1:niter
i = 1; % Testing...
    %%% Play Sounds
    % Make Sounds 
    newstatus(handles,sprintf('%d: Playing Sound',i));
    [stim1, stim2] = golay(order);
    %stim1 = stim1 * amplitude;      % amplitude scaling
    %stim2 = stim2 * amplitude;

    % Calc values
    if i == 1 % should be same on all iters, don't my waste time
        stimpts = length(stim1);
        totalpts = stimpts + length(stim2);
        dur_1 = stimpts/fs;
        dur   = totalpts/fs;
        
        % get list of good freqs
        Ny = fs/2;
        freq = 0 : Ny/((stimpts/2)-1) : Ny;
        lbin = find(freq>=1000,1);
        hbin = find(freq<=70000,1,'last');
    end

    % Play sounds
for j=1:niter
    [rec1(:,j),pbstart1,absrecpos1,recstart1] = playsound(handles,stim1,dur_1);
    pause(0.2)
    [rec2(:,j),pbstart2,absrecpos2,recstart2] = playsound(handles,stim2,dur_1);
    pause(0.5)
end
rec1 = mean(rec1,2);
rec2 = mean(rec2,2);
    %rec1 = rec1.*2;
    %rec2 = rec2.*2;
    
     %Plot in time domain
%     axes(handles.time)
%     cla
%     hold on
%     scatter(1:stimpts,stim1,1,'b','filled')
%     scatter(1:stimpts,rec1,1,'r','filled')
%     line([1:stimpts;1:stimpts],[stim1;rec1],'Color','r')
%     xlim([1,stimpts])
%     ylim([-1,1])
%     hold off
    
    % Calc impulse functions
    newstatus(handles,sprintf('%d: Calculating Impulse Response',i));
    [Imp1,TF1,frequs] = gimpulse(stim1,stim2,rec1,rec2);
    % eliminate negative freqs and handle extreme freqs
    %TF1 = TF1(1:stimpts/2);
    %TF1(1:lbin) = ones(size(lbin))/10000;
    %TF1(hbin:end)=ones(1,length(TF1)-hbin+1)/10000;

    % Plot in freq domain
    %axes(handles.frequency)
    %plot(freq,20*log10(abs(TF1)));
    %xlim([lowfreq,highfreq])
    
        % Plot impulse responses
    axes(handles.time)
    cla
    scatter(([1:length(Imp1)])/fs,Imp1);
    
    axes(handles.frequency)
    semilogx(frequs,TF1)
    xlabel('Frequency [Hz]')
    ylabel('Magnitude [dB]')
    
    
    %Test
    %kernelsize = length(TF1);
    filter_size = 256;
    
    % lag by 1 so is causal?
    
    % Calculate Filter
    newstatus(handles,sprintf('%d: Calculating Filter',i));
    %lsq_gain = [-80, rolloffamp, 0, 0, rolloffamp, -80];
    lsq_freq = [lowrolloff, lowfreq, highfreq, highrolloff];
    filt_TF1 = lsqinv3(TF1',kernelsize, lsq_freq, fs, 0, lsq_gain)';
    [filt_1,MSE1] = lsq2(rec1,stim1,filter_size);
    [filt_2,MSE2] = lsq2(rec2,stim2,filter_size);
    % Plot product of 
        
    %axes(handles.frequency)
    freq_filt = 0:Ny/(kernelsize-1):Ny;
    FT1 = fft(TF1,nfft);
    FTinv1 = fft(filt_TF1,kernelsize * 2);
        
    axes(handles.frequency)
    cla
    hold on
    plot(freq_filt,20*log10(abs(FT1(1:kernelsize))),'r')
    plot(freq_filt,20*log10(abs(FTinv1(1:kernelsize))),'b')
    plot(freq_filt,20*log10(abs(FTinv1(1:kernelsize)))+ 20*log10(abs(FT1(1:kernelsize))),'g')
    xlim([lowfreq,highfreq])
    hold off

%end


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

function [a,b] = golay(n)
%
% create a pair of golay codes of order n
% use assorted randomizations which preserve 
% complementary series properties
%
%IN: n = order of golay sequences
%
%OUT: a = golay sequence 1
%     b = golay sequence 2
%
%note : xcorr(a) + xcorr(b) = delta_fn
%
% Anthony Leonardo, Caltech, 3/01
%

a1 = sign(unifrnd(-1,1,1,1));
a2 = sign(unifrnd(-1,1,1,1));

a = [a1 a2];
b = [a1 -a2];

for i = 2:n
	a2 = [a b];
	b = sign(unifrnd(-1,1,1,1))*[a -b];
	a = sign(unifrnd(-1,1,1,1))*a2;
	if sign(unifrnd(-1,1,1,1)) > 0
		a = fliplr(a);
	end;
	if sign(unifrnd(-1,1,1,1)) > 0
		b = fliplr(b);
	end;
end;

function [I,H,frequs] = gimpulse(a,b,aout,bout)
%
% get the golay-probe based impulse response
%
%IN: a = input to system a
%    b = input to system b
%    aout = output from system a
%    bout = output from system b
%
%OUT: I = impulse response
%     H = transfer function
%
% Anthony Leonardo, Caltech, 3/01
%
% [I,H] = gimpulse(a,b,aout,bout)
%
global pref
 
L = length(a);
M = length(aout);

% de-mean data
a = a - mean(a);
b = b - mean(b);
aout = aout - mean(aout);
bout = bout - mean(bout);

% Compute impulse response
% https://ccrma.stanford.edu/realsimple/imp_meas/imp_meas.pdf
I = fftfilt(a(L:-1:1),aout) + fftfilt(b(L:-1:1),bout);
I = I / (2*L);

% Last ten samples are always fucked up and I don't know why
I(end-9:end) = 0;

% Rescale
I = I/max(abs(I));

% Make freq mag response

frequs = linspace(0,pref.fs/2,L/2+1);
H = fft(I);
H = 20*log10(abs(H(1:L/2+1)));


function b=lsqinv3(h,L,Wn,Fs,pflag,gain)
% LSQINV3 Least square inverse filter design with low pass or band pass filtering
%        b=LSQINV3(h,L,Wn,Fs,pflag, gain)
%        input parameters:
%          h      --> measured impulse response
%          L      --> length of desired impulse response
%	   Wn     --> corner frequencies of bp filter in Hz
%		[lower_stop_band, lower_pass_band, upper_pass_band, upper_stop_band]
%
%			default: no filtering
%	   Fs     --> sampling frequency in Hz (Default Fs = 30000)
%	   pflag     --> plot of transfer functions and impulse responses
%		(1: yes, 2: no, default no)
%
%     gain dB attenuation at each band of Wn:
%		[lower_stop_band, lower_pass_band, passband,passband, upper_pass_band, upper_stop_band]
%
%        output parameters:
%          b --> impulse response of inverse filter
%
% LSQINV solves the set of linear equations
%
%  |rhh(0) rhh(1)   ... rhh(N)  ||b(0)|   |rdh(0)  |
%  |rhh(1) rhh(0)   ... rhh(N-1)||b(1)|   |rdh(1)|
%  | .      .            .      || .   | = | .    |
%  | .      .            .      || .   |   | .    |
%  |rhh(N) rhh(N-1) ... rhh(0)  ||b(N)|   |rdh(N)|
%
% where h(k)=0 for k<0 and rhh(k) is the autocorrelation function of h

% References:   Proakis & Manolakis
%               Digital Signal Processing
%               2nd ed.
%
% low pass filtering response bandpass filtering
% before linear equations are solved
% solving rdh(n) instead of h(1)

L_h = length(h);

d = zeros(1,L);
h = [h; zeros(L-L_h,1)];
b = zeros(L,1);



  f = [ 0, Wn(1), Wn(2), Wn(3), Wn(4), Fs/2];             % corner frequs
  m_dB = gain;
  m = 10.^(m_dB/20);
  m(6) = 0;
  d = firls(L-1, 2*f/Fs, m);
  %d = firpm(17,2*f/Fs,m);
  [h,w] = freqz(d,1,L);
  
  %d = d(1:L);      % added 3/25/03 because firls increases order by one if gain at Nyquist ~=0

  d = d(:);         



rhh=xcorr(h);
rhh=rhh(L_h:L_h+L-1);
R = toeplitz(rhh);
m = size(d,1);
%rdh = xcorr(d,flipud(h(:,nc)));
rdh = xcorr(d,h);
rdh = rdh(m:end);

b = R\rdh;
 


    
function [B,MSE] = lsq2(f,d,M)
% Calculate least-squares inverse filter
% From https://ocw.mit.edu/courses/mechanical-engineering/2-161-signal-processing-continuous-and-discrete-fall-2008/study-materials/lsqfit.pdf
% f - observed samples
% d - desired samples
% M - filter order
% B - filter coefficients
% MSE - mean-square-error

N = length(f);

phiff = xcorr(f);
phifd = xcorr(d,f);

rff = phiff(N:N+M-1);
R   = toeplitz(rff);
P=phifd(N:N+M-1);

B=inv(R)*P';
phidd=xcorr(d);
MSE=phidd(N)-P*B;


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
buffSize        = 2048; % Need big buffer (golay order 11) for some reason
buffSize = [];
reqlatencyclass = 3;      % Try for lowest latency
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

function [recorded,pbstart,absrecpos,recstart] = playsound(handles,sound,dur)
% Dur in s
% Have to pad with silence due to latency of soundcard, we
% recover the real recording later.
sound_pad = [zeros(2000,1)',sound];
pahandle=handles.audio;
PsychPortAudio('FillBuffer', pahandle, sound_pad); % fill playback buffer
PsychPortAudio('GetAudioData', pahandle, dur*5); % Preallocate a recording buffer

nreps=2;when=0;waitForStart=1;
pbstart = PsychPortAudio('Start', pahandle,nreps,when,waitForStart);
waitForEndOfPlayback=1; %1: Wait until playback finished to stop
PsychPortAudio('Stop', pahandle, waitForEndOfPlayback);

% Retrieve audio data
[recorded, absrecpos, overflow, recstart] = PsychPortAudio('GetAudioData', pahandle);
if overflow>0
   newstatus(handles, 'Warning! Overflow in recording!')
end

% Find actual start of recording
% mic is rl quiet before actual start, no sample exceeds 0.01, but first
% sample of recording usually does.
% We use diff because the mic can sometimes float above/below zero, a
% strong departure is a better indicator than an absolute value
startind = find(abs(diff(recorded))>0.008,1)+2; % Want one sample after so filter is causal
recorded = recorded(startind:startind+length(sound)-1);

% Normalize recorded sound like time sound
%recorded = ((recorded-min(recorded))./(max(recorded)-min(recorded)))-0.5;
%recorded = recorded-mean(recorded(:));
%recorded = recorded/(std(recorded(:))*2);

function newstatus(handles,string)
cur_strings = cellstr(handles.status_list.String);
new_strings = [cur_strings;{string}];
set(handles.status_list,'String',new_strings)
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% Useless shit that no one cares about
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function highfreq_val_Callback(hObject, eventdata, handles)

function highfreq_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to highfreq_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lowfreq_val_Callback(hObject, eventdata, handles)
% hObject    handle to lowfreq_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%maxdur = str2num(get(hObject,'String'));
%handles.time.XLim(2) = maxdur;
%handles.measured.XLim(2)    = maxdur;
%pahandle = handles.audio;
%PsychPortAudio('GetAudioData', pahandle, maxdur/1000);

function lowfreq_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowfreq_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function highrolloff_val_Callback(hObject, eventdata, handles)

function highrolloff_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to highrolloff_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lowrolloff_val_Callback(hObject, eventdata, handles)
% hObject    handle to lowrolloff_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

maxfreq = str2num(get(hObject,'String'));
handles.time.YLim(2) = maxfreq;
handles.measured.YLim(2)    = maxfreq;
handles.frequency.XLim(2)      = maxfreq;

function lowrolloff_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowrolloff_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function order_val_Callback(hObject, eventdata, handles)
% hObject    handle to order_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%minfreq = str2num(get(hObject,'String'));
%handles.time.YLim(1) = minfreq;
%handles.measured.YLim(1)    = minfreq;
%handles.frequency.XLim(1)      = minfreq;

function order_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to order_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rolloffamp_val_Callback(hObject, eventdata, handles)

function rolloffamp_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rolloffamp_val (see GCBO)
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

function kernel_val_Callback(hObject, eventdata, handles)

function kernel_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kernel_val (see GCBO)
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


% --- Executes when user attempts to close complexcalibration_golay.
function complexcalibration_golay_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to complexcalibration_golay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
PsychPortAudio('Stop', handles.audio);
PsychPortAudio('Close', handles.audio);
end
% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on selection change in status_list.
function status_list_Callback(hObject, eventdata, handles)
% hObject    handle to status_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns status_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from status_list


% --- Executes during object creation, after setting all properties.
function status_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to status_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function complexcalibration_golay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to complexcalibration_golay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
