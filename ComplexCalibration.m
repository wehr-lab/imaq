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

% Last Modified by GUIDE v2.5 11-Jan-2017 18:10:27

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

% UIWAIT makes ComplexCalibration wait for user response (see UIRESUME)
% uiwait(handles.figure1);


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
minfreq  = str2num(handles.minfreq_val.String);
maxfreq  = str2num(handles.maxfreq_val.String);
dur      = str2num(handles.dur_val.String);
am_vel   = str2num(handles.am_val.String);
fm_dens  = str2num(handles.fm_val.String);
depth    = str2num(handles.mod_val.String);
niter    = str2num(handles.niter_val.String);
converge = str2num(handles.converge_val.String);

% just testing...
fs = 192000;

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



%for i = 1:niter
   handles.status_text.String = 'Generating Sound';
   drawnow
   
   snd = genripple(am_vel,fm_dens,depth,dur,fs,freq);
   
   handles.status_text.String = 'Generating Spectrogram';
   axes(handles.synthesized);
   spectrogram(snd,winsize,noverlap,nfft,fs,'yaxis');
   %set(gca,'YScale','log');
   ylim([minfreq/1000,maxfreq/1000]);
   colorbar('off');
   
   handles.status_text.String = 'Generating Periodogram';
   axes(handles.filter);
   periodogram(snd,[],nfft,fs);
   xlim([minfreq/1000,maxfreq/1000]);
   
   handles.status_text.String = '';
   drawnow
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
