function varargout = StandaloneCalibrateSpeaker(varargin)

%Help for StandaloneCalibrateSpeaker
%
% Forked to imaq from djmaus by JLS, 1/9/2017
% Minor compatibility changes.
%
%
% Description:
% This function creates an equalization table of attenuation values for
% amplitude/frequency pairs in a file called "calibration.mat" in the
% djmaus root directory.
% This table is esentially the inverse of the speaker's frequency response
% curve. This look-up table is then used by djmaus to equalize sounds
% such tones and white noise. This process is called equalization.
%
% -------------------------------------------------------------------------
%
% Directions for free-field calibration:
% First hook up the output from the B&K microphone (BNC) to the soundcard IN1 (XLR)
% using a BNC-XLR adaptor.  Then turn down the amplifier (e.g. Sampson) that
% the speaker is connected to in order to avoid damaging
% the speaker.  The value in the target dB field will be used to calculate
% the needed attenuation.  Usually you need to loop the process a number of times
% to produce a more even attenuation curve.  To change the number of loops put
% desired number in "num Loops" field.  The value in the maxfreq field
% should be anywhere from 32000 to 80000, depending on the speakers ability
% to play high frequency tones.
% After an initial calibration, turn up the amplifier a
% small amount if target dB value is not being achieved.
%
% -freqs per octave: frequency resolution to use between the minfreq and maxfreq
%   values. Note that an extra frequency is added to calibrate whitenoise.
% -save: whether to overwrite the saved calibration table. Unclick save
%   (the default) to just check the current calibration status. If
%   you want to re-calibrate (and overwrite the saved calibration!) then
%   click save on.
%   Note: you have to Run at least once with save on before it will save to a file
% -num loops: how many times to iterate.
% -convergence: inverse filter convergence factor. Controls how fast to converge
%   towards inverse filter (defaults to 1.0, reduce to maybe .9 or lower if oscillating
%   around target amplitude)
% -reset_atten: how much to reset the minimum attenuation towards zero each
%   iteration. Start with 0. If min attenuation starts to
%   accumulate (e.g. more than 5 dB) set it to 0.5 or so. Setting to 1
%   will work faster but can sometimes oscillate. The min atten is reported
%   in the message window at completion along with other data.
% -Reset: click this button to erase data stored in the calibration file
%   and start over from scratch.
% -mic_sensitivity: the gain setting on the B&K Nexus amplifier. Defaults to
%   1 V/Pa. If for some reason you change the B&K gain (e.g. to .316 V/Pa)
%   then set this value accordingly.
% -override_atten: do not use this, leave it set to zero. Advanced users:
%   this overrides atten by directly subtracting the specified value from atten.
% Notes:
% -The soundcard DeviceID should be set automatically to the ASIO device taken from djPrefs,
% but you can override it if you want to. Input Ch and Output Ch should be
% left set to 1.
% -Negative attenuation values are set to 0.
% -Sampling rate and number of channels are set for you, so you don't need to worry
%   about it
% -If you change the min, max, or number of frequencies, the table will be
%   erased and calibration will start over from scratch.
% - pure tone amplitude is computed from the peak in the power spectrum.
%   white noise amplitude is computed as rms. White noise data is high-pass
%   filtered at 500 hz to remove low-frequency ambient noise contamination
%
% -------------------------------------------------------------------------

global user
global pref
if ~nargin
    [ok, user]=Login;
    if ~ok
        return
    end
end

warning off MATLAB:Axes:NegativeDataInLogAxis
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @StandaloneCalibrateSpeaker_OpeningFcn, ...
    'gui_OutputFcn',  @StandaloneCalibrateSpeaker_OutputFcn, ...
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

% --- Executes just before StandaloneCalibrateSpeaker is made visible.
function StandaloneCalibrateSpeaker_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to StandaloneCalibrateSpeaker (see VARARGIN)

global pref;
global user;
Prefs(user);

% Choose default command line output for StandaloneCalibrateSpeaker
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using StandaloneCalibrateSpeaker.
if strcmp(get(hObject,'Visible'),'off')
    t=1:1000;t=t/1000;
    plot(t,sin(2*pi*10*t));
end

%Initialize default parameters, daq, etc.
InitParams(handles)
InitSound(handles)

%why is this done yet again, here? commenting out 090216
% DeviceID = 1+get(handles.DeviceID, 'Value');
% userdata=get(handles.figure1, 'userdata');
% userdata.DeviceID=DeviceID;
% set(handles.figure1, 'userdata', userdata);
% InitializeInputCh(handles)
% InitSoundIn(handles)
%
% DeviceID = get(handles.DeviceID, 'Value')-1;
% userdata=get(handles.figure1, 'userdata');
% userdata.DeviceID=DeviceID;
% set(handles.figure1, 'userdata', userdata);
% InitializeOutputCh(handles)
% InitSoundOut(handles)

% UIWAIT makes StandaloneCalibrateSpeaker wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = StandaloneCalibrateSpeaker_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pref
axes(handles.axes1);
cla;
axes(handles.axes2);
cla;
userdata=get(handles.figure1, 'userdata');
samplingrate=userdata.samplingrate;
atten=userdata.atten;
logspacedfreqs=userdata.logspacedfreqs;
numfreqs=length(logspacedfreqs);
logspacedfreqs=[-1 logspacedfreqs]; %add white noise at end
maxfreq=max(logspacedfreqs);
target_amp=str2num(get(handles.target_amp, 'string'));
convergence=str2num(get(handles.convergence, 'string'));
override_atten=str2num(get(handles.override_atten, 'string'));
reset_atten=str2num(get(handles.reset_atten, 'string'));
save1=get(handles.save,'value');
h=handles.minfreq;
minfreq=str2num(get(h, 'string'));
h=handles.freqsperoctave;
freqsperoctave=str2num(get(h, 'string'));

Message('Running...', handles)
set(hObject, 'background', 'r', 'String', 'Running...', 'enable', 'on')

wb=waitbar(0, 'calibrating...');
set(wb, 'pos', [870 700 270  50])
num_loops=str2num(get(handles.numLoops, 'string'));

for lp = 1:num_loops
    sampleLength=500;   % sample length in ms
    tonedur=sampleLength+300;
    %                 numfreqs=userdata;
    %                 maxfreq=GetParam(me, 'maxfreq');
    %                 minfreq=GetParam(me, 'minfreq');
    %                 if minfreq==-1; error('minfreq not supposed to be -1'); end
    %                 logspacedfreqs = logspace( log10(minfreq) , log10(maxfreq) , numfreqs );
    if isempty (atten)
        atten=zeros(size(logspacedfreqs));
    elseif length(atten)~=length(logspacedfreqs)
        atten=zeros(size(logspacedfreqs));
    end
    
    for i=1:length(logspacedfreqs) %play tones
        running=get(hObject, 'value');
        if ~running %user pressed run again
            break
        end
        
        Message(sprintf('playing tone %d/%d (%.1f Hz)',i, numfreqs+1, logspacedfreqs(i)), handles)
        tonefreq=logspacedfreqs(i);
        %start playing tone
        audiodata=PlayCalibrationTone(tonefreq, tonedur, handles);
        BKsensitivity=str2num(get(handles.BKsensitivity, 'string'));
        ScaledData=detrend(audiodata(1,:), 'constant'); %in volts
        %high pass filter a little bit to remove rumble
        %for display purposes only
        [b,a]=butter(1,100/(samplingrate), 'high');
        % Display trace.
        ax1=handles.axes1;
        ax2=handles.axes2;
        t=1:length(ScaledData);
        t=1000*t/samplingrate;
        plot(ax1,t,filtfilt(b,a,ScaledData));
        xlabel(ax1, 'time (ms)');
        ylabel(ax1, 'Microphone Voltage (V)');
        
        %estimate frequency
        hold(ax2, 'on')
        xlim(ax2, [0 1.25*maxfreq])
        [Pxx,f] = pwelch(ScaledData,[],[],[],samplingrate);
        fmaxindex=(find(Pxx==max(Pxx(10:end)))); %skip freqs<210hz (for 192kHz Fs)
        fmaxindex=fmaxindex(1);
        fmax=round(f(fmaxindex));
        Message(sprintf('est. freq: %d', fmax), handles);
        
        c=repmat('rgbkm', 1, ceil(numfreqs/5)+1);
        semilogy(ax2, f(10:end), Pxx(10:end), c(i))
        semilogy(ax2, f(fmaxindex), Pxx(fmaxindex), ['o',c(i)])
        xlabel(ax2, 'Frequency, kHz');
        xt=get(ax2, 'xtick');
        set(ax2, 'xticklabel', round(10*xt/1000)/10)
        
        %if GetXonarDevice & isempty(GetAsioLynxDevice)
        %    fudgefactorTone=-1.07;
        %    fudgefactorWN=0.79;
        %if GetAsioLynxDevice  & isempty(GetXonarDevice)
            fudgefactorTone=18.7;
            fudgefactorWN=9.3;
        %else
        %    error ('what soundcard are we using for recording?')
        %end
        ylabel(ax2, 'PSD');
        if tonefreq==-1
            %estimate amplitude -- RMS method
            % high-pass filtering at 500 hz (mw 01-30-09)
            %ai2SampleRate
            [b,a]=butter(5, 500/(samplingrate/2), 'high');
            Vrms=sqrt(mean(filtfilt(b,a,ScaledData).^2));
            db=dBSPL(Vrms, BKsensitivity);
            %fudge factor to achieve 94dB for B&K calibrator:
            db=db+fudgefactorWN;
            Message( sprintf('estimated noise amp: %.2f dB', db), handles);
        else
            %estimate amplitude -- Pxx method
            %                 fidx=closest(f, tonefreq);
            db=dBPSD(Pxx(fmaxindex), BKsensitivity);%should return 94 for B&K calibrator
            %fudge factor to achieve 94dB for B&K calibrator:
            db=db+fudgefactorTone;
            %db=dBPSD(Pxx(fidx), GetParam(me, 'mic_sensitivity'));
            Message( sprintf('estimated tone amp: %.2f dB', db), handles);
            %                             pause(.35)
        end
        FMAX(i)=fmax;
        DB(i)=db;
        waitbar(((lp-1)*numfreqs+i)/(num_loops*numfreqs), wb)
    end %play tones
    if ~running %user pressed run again
        break
    end
    
    Message('measured freqs: ', handles)
    Message(sprintf('%d ',FMAX), handles)
    Message('measured dB: ', handles)
    Message(sprintf('%.1f ',DB), handles)
    
    
    %since WN shows up as -1, which doesn't work with log plots, sticking it
    %just below minfreq
    logspacedfreqs_forplotting=logspacedfreqs;
    if logspacedfreqs_forplotting(1)==-1
        logspacedfreqs_forplotting(1)=logspacedfreqs_forplotting(2)^2/logspacedfreqs_forplotting(3);
        semilogx(ax1, logspacedfreqs_forplotting(1), DB(1), 'r*')
        text(ax1, logspacedfreqs_forplotting(1), DB(1), 'WN')
        hold(ax1, 'on')
    end
    semilogx(ax1, logspacedfreqs_forplotting, DB, '-o')
    ylim(ax1, [.9*min(DB) 1.1*max(DB)])
    hold(ax1, 'off')
    xlabel(ax1, 'Frequency, kHz');
    ylabel(ax1, 'dB SPL');
    xlim(ax1, [0 1.25*maxfreq])
    grid(ax1, 'on')
    set(ax1, 'xtick', logspacedfreqs)
    set(ax1, 'xticklabelmode', 'auto')
    xt=get(ax1, 'xtick');
    set(ax1, 'xticklabel', round(10*xt/1000)/10)
    
    %  create inverse filter
    %atten=DB-min(DB); %this is the inverse filter
    atten=convergence*(DB-target_amp); %this is the inverse filter
    %atten(atten<0)=0;
    
    %apply over ride value if any
    atten=atten-override_atten;
    
    % iteratively store calibration data
    if save1
        stored_atten=userdata.atten;
        if length(atten)==length(stored_atten)
            atten=atten+stored_atten; %iteratively add to stored calibration
        end
        reset_amount=reset_atten;
        atten=atten-reset_amount*min(atten); %reset min atten towards 0 to avoid saturating
        atten(atten<0)=0;
        atten % print for debugging
        cd(pref.calibration)
%         try
%             cal=load('calibration.mat');
%         catch % #ok
%             cal=[];
%         end
        timestampstr=['last saved ', datestr(now)];
        save calibration logspacedfreqs timestampstr DB  atten minfreq maxfreq freqsperoctave
        userdata.atten=atten;
        userdata.DB=DB;
        set(handles.figure1, 'userdata', userdata);
        Message('saved calibration data to file', handles)
        
    end
    Message([sprintf('\nmin dB %.2f', min(DB)), sprintf('\nmax dB %.2f', max(DB)), ...
        sprintf('\nmean dB %.2f', mean(DB)), sprintf('\nstd dB %.2f', std(DB)), ...
        sprintf('\nmin atten %.2f', min(atten)),  ...
        sprintf('\nmax atten %.2f', max(atten)), ...
        sprintf('\npref.maxSPL is %d',pref.maxSPL )], handles)
    StdDB(lp)=std(DB);
end %loop_num

%stop and turn off Run button
set(hObject, 'background', [0 1 0], 'String', 'Run','enable', 'on', 'value', 0)

close(wb)
Message('Done', handles)




% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
    ['Close ' get(handles.figure1,'Name') '...'],...
    'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});


% --- Executes on button press in Reset.
function Reset_Callback(hObject, eventdata, handles)
% hObject    handle to Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Reset(handles)

function Reset(handles)
userdata=get(handles.figure1, 'userdata');
userdata.atten=[];
userdata.DB=[];
set(handles.figure1, 'userdata', userdata);
% InitParams(handles)
cla(handles.axes1)
cla(handles.axes2)

function minfreq_Callback(hObject, eventdata, handles)
% hObject    handle to minfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minfreq as text
%        str2double(get(hObject,'String')) returns contents of minfreq as a double
calculate_logspacedfreqs(handles);
Reset(handles)

% --- Executes during object creation, after setting all properties.
function minfreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxfreq_Callback(hObject, eventdata, handles)
% hObject    handle to maxfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxfreq as text
%        str2double(get(hObject,'String')) returns contents of maxfreq as a double
calculate_logspacedfreqs(handles);
Reset(handles)

% --- Executes during object creation, after setting all properties.
function maxfreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function freqsperoctave_Callback(hObject, eventdata, handles)
% hObject    handle to freqsperoctave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freqsperoctave as text
%        str2double(get(hObject,'String')) returns contents of freqsperoctave as a double
calculate_logspacedfreqs(handles);
Reset(handles)

% --- Executes during object creation, after setting all properties.
function freqsperoctave_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freqsperoctave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function target_amp_Callback(hObject, eventdata, handles)
% hObject    handle to target_amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of target_amp as text
%        str2double(get(hObject,'String')) returns contents of target_amp as a double
Reset(handles)

% --- Executes during object creation, after setting all properties.
function target_amp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to target_amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save



function convergence_Callback(hObject, eventdata, handles)
% hObject    handle to convergence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of convergence as text
%        str2double(get(hObject,'String')) returns contents of convergence as a double


% --- Executes during object creation, after setting all properties.
function convergence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to convergence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function reset_atten_Callback(hObject, eventdata, handles)
% hObject    handle to reset_atten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reset_atten as text
%        str2double(get(hObject,'String')) returns contents of reset_atten as a double


% --- Executes during object creation, after setting all properties.
function reset_atten_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reset_atten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function override_atten_Callback(hObject, eventdata, handles)
% hObject    handle to override_atten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of override_atten as text
%        str2double(get(hObject,'String')) returns contents of override_atten as a double


% --- Executes during object creation, after setting all properties.
function override_atten_CreateFcn(hObject, eventdata, handles)
% hObject    handle to override_atten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

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

function InitSound(handles)
global pref
% InitializePsychSound(0); In imaq this is done by Prefs
%InitializeInputCh(handles) This gets done in InitParams
%InitializeOutputCh(handles)
userdata=get(handles.figure1, 'userdata');
SoundFs=userdata.samplingrate;
DeviceID=userdata.DeviceID;
outputCh=userdata.OutputCh;
reqlatencyclass =1; %Level 1 (the default) means: Try to get the lowest latency that is possible under the constraint of reliable playback, freedom of choice for all parameters and interoperability with other applications
numChan=[1,1]; % Two output (L/R, 1 input)
buffSize=[]; %use default
Mode=3; %full duplex: simultaneous capture and playback

%hack mw 10.13.2016 to use other soundcard for input
%Mode=1;

%stop and close
try
    PsychPortAudio('Stop', handles.pahandle);
    PsychPortAudio('Close', handles.pahandle);
end

try pahandle = PsychPortAudio('Open', DeviceID, Mode, reqlatencyclass, SoundFs, numChan, buffSize);
 %hack mw 10.13.2016:
 %paInhandle = PsychPortAudio('Open', 0, 2, reqlatencyclass, SoundFs, numChan, buffSize);
    %runMode = 0; %default, turns off soundcard after playback
    runMode = 1; %leaves soundcard on (hot), uses more resources but may solve dropouts? mw 08.25.09: so far so good.
    PsychPortAudio('RunMode', pahandle, runMode);
    % Preallocate an internal audio recording  buffer with a capacity of 10 seconds:
    PsychPortAudio('GetAudioData', pahandle, 10);%hack mw 10.13.2016
    
    userdata.pahandle=pahandle;
    %userdata.paInhandle=paInhandle;%hack mw 10.13.2016
    set(handles.figure1, 'userdata', userdata);
    Message('Initialized Sound', handles)
    
catch
    Message(sprintf('Error: could not open soundcard Device'), handles);
    set(handles.Run, 'enable', 'off')
end



function InitParams(handles)
global pref
minfreq=1000;
maxfreq=2e3;
freqsperoctave=4;
target_amp=70;
save=0;
convergence=1;
reset_atten=0;
override_atten=0;


try %try to load existing calibration
    cal=load(fullfile(pref.calibration, 'calibration.mat'));
    atten=cal.atten;
    DB=cal.DB;
    logspacedfreqs=cal.logspacedfreqs;
    semilogx(handles.axes1, logspacedfreqs, DB, '-o')
    ylabel('stored dB')
    xlabel('frequency, Hz')
    grid on
    minfreq=cal.minfreq;
    maxfreq=cal.maxfreq;
    freqsperoctave=cal.freqsperoctave;
    Message('succesfully loaded previously saved calibration', handles);
catch
    Message('failed to load previously saved calibration', handles);
    atten=[];
end

numoctaves=log2(maxfreq/minfreq);
logspacedfreqs=minfreq*2.^([0:(1/freqsperoctave):numoctaves]);
newmaxfreq=logspacedfreqs(end);
numfreqs=length(logspacedfreqs);
if maxfreq~=newmaxfreq
    Message(sprintf('note: could not divide %d-%d Hz evenly into exactly %d frequencies per octave', minfreq, maxfreq, freqsperoctave), handles)
    Message(sprintf('using new maxfreq of %d to achieve exactly %d frequencies per octave', round(newmaxfreq), freqsperoctave), handles)
    maxfreq=newmaxfreq;
end

h=handles.minfreq;
set(h, 'string', minfreq)
h=handles.maxfreq;
set(h, 'string', maxfreq)
h=handles.freqsperoctave;
set(h, 'string', freqsperoctave)
h=handles.target_amp;
set(h, 'string', target_amp)
h=handles.save;
set(h, 'value', save)
h=handles.convergence;
set(h, 'string', convergence)
h=handles.reset_atten;
set(h, 'string', reset_atten)
h=handles.override_atten;
set(h, 'string', override_atten)


userdata=get(handles.figure1, 'userdata');
userdata.samplingrate= pref.fs;
userdata.DeviceID= pref.dev_id;

set(handles.figure1, 'userdata', userdata);
InitializeInputCh(handles)
InitializeOutputCh(handles)
userdata.InputCh=handles.InputCh;
userdata.OutputCh=handles.OutputCh;

userdata.atten=atten;
userdata.DB=[];
userdata.logspacedfreqs=logspacedfreqs;
set(handles.figure1, 'userdata', userdata);



function InitializeInputCh(handles)
userdata=get(handles.figure1, 'userdata');
devs = PsychPortAudio('GetDevices'); %devs is zero-indexed, hence the confusion
DeviceID = userdata.DeviceID;
numchannels=devs(DeviceID+1).NrInputChannels;
for i = 1:numchannels
    ChString{i}=sprintf('%d', i);
end
if numchannels==0 ChString='No Input Ch!!';
    Message('Error! Bad Input Device ID', handles)
end
set(handles.InputCh, 'String', ChString);
set(handles.DeviceID, 'Value', DeviceID+1)

function InitializeOutputCh(handles)
userdata=get(handles.figure1, 'userdata');
devs = PsychPortAudio('GetDevices');
DeviceID = userdata.DeviceID;
numchannels=devs(DeviceID+1).NrOutputChannels;
for i = 1:numchannels
    ChString{i}=sprintf('%d', i);
end
if numchannels==0
    ChString='No Output Ch!!';
    Message('Error! Bad Output Device ID', handles)
end
set(handles.OutputCh, 'String', ChString);
set(handles.DeviceID, 'Value', DeviceID+1)


% --- Executes on selection change in InputCh.
function InputCh_Callback(hObject, eventdata, handles)
% hObject    handle to InputCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns InputCh contents as cell array
%        contents{get(hObject,'Value')} returns selected item from InputCh
InputCh = get(handles.InputCh, 'Value');
userdata=get(handles.figure1, 'userdata');
userdata.InputCh=InputCh;
set(handles.figure1, 'userdata', userdata);
InitSound(handles)

% --- Executes during object creation, after setting all properties.
function InputCh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InputCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in OutputCh.
function OutputCh_Callback(hObject, eventdata, handles)
% hObject    handle to OutputCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns OutputCh contents as cell array
%        contents{get(hObject,'Value')} returns selected item from OutputCh
OutputCh = get(handles.OutputCh, 'Value');
userdata=get(handles.figure1, 'userdata');
userdata.OutputCh=OutputCh;
set(handles.figure1, 'userdata', userdata);
InitSound(handles)

% --- Executes during object creation, after setting all properties.
function OutputCh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OutputCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in DeviceID.
function DeviceID_Callback(hObject, eventdata, handles)
% hObject    handle to DeviceID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns DeviceID contents as cell array
%        contents{get(hObject,'Value')} returns selected item from DeviceID
DeviceID = get(handles.DeviceID, 'Value')-1;
userdata=get(handles.figure1, 'userdata');
userdata.DeviceID=DeviceID;
set(handles.figure1, 'userdata', userdata);
InitializeInputCh(handles)
InitializeOutputCh(handles)
InitSound(handles)


% --- Executes during object creation, after setting all properties.
function DeviceID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DeviceID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
InitializePsychSound(0)
devs = PsychPortAudio('GetDevices');
for i = 1:length(devs)
    deviceString{i}=sprintf('%d: %s: %s', devs(i).DeviceIndex, devs(i).HostAudioAPIName, devs(i).DeviceName);
end
set(hObject, 'String', deviceString);
set(hObject, 'Value',7); %default Output Device ID, note devs is 0-indexed



function Message_Callback(hObject, eventdata, handles)
% hObject    handle to Message (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Message as text
%        str2double(get(hObject,'String')) returns contents of Message as a double


% --- Executes during object creation, after setting all properties.
function Message_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Message (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    % % Stop capture:
    %PsychPortAudio('Stop', handles.paInhandle);
    % % Close the audio device:
    %PsychPortAudio('Close', handles.paInhandle);
    % % Stop playback:
end
try
    PsychPortAudio('Stop', handles.pahandle);
    % % Close the audio device:
    PsychPortAudio('Close', handles.pahandle);
end
PsychPortAudio('Close') %If pahandle is omitted, all audio devices will be closed and the driver will shut down.

function audiodata=PlayCalibrationTone(tonefreq, tonedur, handles);
userdata=get(handles.figure1, 'userdata');
target_amp=str2num(get(handles.target_amp, 'string'));
param.frequency=tonefreq; %hz
param.amplitude=target_amp;
nstimchans=1;

logspacedfreqs=userdata.logspacedfreqs;
atten=userdata.atten;
if tonefreq==-1 findex=1;
else
    findex=find(logspacedfreqs<=tonefreq, 1, 'last')+1;
end
if isempty(atten)
    attenuation=0;
else
    attenuation=atten(findex);
end
param.amplitude=param.amplitude-attenuation;

param.duration=tonedur; %ms
param.ramp=10;
samplerate=userdata.samplingrate;
if tonefreq==-1
    samples=MakeWhiteNoise(param, samplerate);
else
    samples=MakeTone(param, samplerate);
end
samples=reshape(samples, nstimchans, length(samples)); %ensure samples are a row vector
%samples(2,:)=0.*samples;

pahandle=userdata.pahandle;
%paInhandle=userdata.paInhandle;
PsychPortAudio('FillBuffer', pahandle, samples); % fill buffer
nreps=1;
when=0; %use this to start immediately
waitForStart=0;
PsychPortAudio('Start', pahandle,nreps,when,waitForStart);
%PsychPortAudio('Start', paInhandle,nreps,when,waitForStart);

% Stop capture:
waitForEndOfPlayback=1; %'waitForEndOfPlayback' - If set to 1, this method will wait until playback of the audio stream finishes by itself.
PsychPortAudio('Stop', pahandle, waitForEndOfPlayback);
%PsychPortAudio('Stop', paInhandle, waitForEndOfPlayback);

% Retrieve pending audio data from the drivers internal ringbuffer:
[audiodata absrecposition overflow cstarttime] = PsychPortAudio('GetAudioData', pahandle);
nrsamples = size(audiodata, 2);
if overflow>0 warning('overflow in RecordTone'); end

function calculate_logspacedfreqs(handles)
h=handles.minfreq;
minfreq=str2num(get(h, 'string'));
h=handles.maxfreq;
maxfreq=str2num(get(h, 'string'));
h=handles.freqsperoctave;
freqsperoctave=str2num(get(h, 'string'));

numoctaves=log2(maxfreq/minfreq);
logspacedfreqs=minfreq*2.^([0:(1/freqsperoctave):numoctaves]);
newmaxfreq=logspacedfreqs(end);
numfreqs=length(logspacedfreqs);
if maxfreq~=newmaxfreq
    Message(sprintf('note: could not divide %d-%d Hz evenly into exactly %d frequencies per octave', minfreq, maxfreq, freqsperoctave), handles)
    Message(sprintf('using new maxfreq of %d to achieve exactly %d frequencies per octave', round(newmaxfreq), freqsperoctave), handles)
    maxfreq=newmaxfreq;
end
h=handles.maxfreq;
set(h, 'string', maxfreq);

userdata=get(handles.figure1, 'userdata');
userdata.logspacedfreqs=logspacedfreqs;
set(handles.figure1, 'userdata', userdata);



function BKsensitivity_Callback(hObject, eventdata, handles)
% hObject    handle to BKsensitivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BKsensitivity as text
%        str2double(get(hObject,'String')) returns contents of BKsensitivity as a double


% --- Executes during object creation, after setting all properties.
function BKsensitivity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BKsensitivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ClearMessages.
function ClearMessages_Callback(hObject, eventdata, handles)
% hObject    handle to ClearMessages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=handles.Message;
set(h, 'String', '');


% --- Executes on button press in Help.
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
help StandaloneCalibrateSpeaker
Message('See command window for help.', handles)

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

function dBSPL=dBSPL(Vrms, sensitivity)
%usage: dBSPL=dBSPL(Vrms, sensitivity)
%inputs: 
%   Vrms: cycle RMS voltage measured from the output of the B&K microphone
%   sensitivity: sensitivity setting on the B&K amplifier, in volts 
%                   (e.g. .1, .316, 1  V/Pa)
%   Vrms = (peak to peak)/sqrt(2)
%output: the power of the sound in dB SPL (re: 20 microPa)

dBSPL=20*log10((Vrms/sensitivity)/20e-6);

function dBPSD=dBPSD(Px, sensitivity)
%usage: dBPSD=dBPSD(Px, sensitivity)
%inputs: 
%   Px: power spectral density of the voltage measured from the output of the B&K microphone
%        this is one element of the Pxx array returned by pwelch
%   sensitivity: sensitivity setting on the B&K amplifier, in volts 
%                   (e.g. .1, .316, 1  V/Pa)
%output: the power of the sound in dB SPL (re: 20 microPa)
%
%algorithm: worked backwards from measurement using the B&K calibrator
%at .316 v/Pa:
%94 dB produces a Pxx of .0076 (.31 Vrms, 1 Vp-p)
%114 dB produces a Pxx of .7620 (3.1 Vrms, 9Vp-p)
%
%at 1 v/Pa:
%94 dB produces a Pxx of .0757 (1 Vrms, 3 Vp-p)
%114 dB produces a Pxx of 5.6720 (9.75 Vrms, 27 Vp-p)
%
%at 3.16 v/Pa:
%94 dB produces a Pxx of .7620 (3.10 Vrms
%114 dB produces a Pxx of 6.96 (clips on scope)

%the Pxx is in units of power(dB)/Hz


switch sensitivity
    case .316
        F=.00757;
    case 1
        F=.0757;
    case 3.16
        F=.757;
    otherwise
        error('expected sensitivity of .316, 1, or 3.16')
end
dBPSD=94+10*log10(Px/F);

