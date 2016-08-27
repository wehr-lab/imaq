function chordtrain=MakeChordTrain(varargin)

global exper pref

% Creates a chord train
% Input
%   nchords         -   number of chords
%   isi             -   inter-stimulus interval, i.e. interval between the
%                       start of previous chord and start of the next chord
%   chordduration   -   duration of an individual chord (ms)
%   frequency       -   frequency components of each chord (Hz)
%   amplitude       -   amplitude of individual chord component !!!!(dB)
%   start           -   start of the first chord after the trigger (ms)
%   duration        -   total duration of the chord train (ms)
%   next            -   inter-train-interval, i.e. when the next
%                       train should follow the previous one (ms)
%   ramp            -   rising/falling edge of individual chords
%   samplerate      -   sampling rate (Hz)
% Output
%   chordtrain      -   requested train
%

chordtrain=[];

if nargin<2
    return;
end

cd(pref.experhome)
cd calibration
cal=load('calibration');

params=varargin{1};
samplerate=varargin{2};
start           =params.start/1000;              % in s
chordduration   =params.chordduration/1000;      % in s
freqs           =params.frequency;
amplitude       =params.amplitude;
nchords         =params.nchords;
isi             =params.isi/1000;                % in s
ramp            =params.ramp;

nfreqs=length(freqs);

train_length=(start+chordduration+(nchords-1)*isi); % in s

%prepare the samples
chordtrain=zeros(ceil(train_length*samplerate),1);

%     sampled_duration=round(chordduration*samplerate);
sampled_start=floor(start*samplerate);
sampled_isi=round(isi*samplerate);

toneparam.duration=chordduration*1000;
toneparam.ramp=0;               % no ramp in individual components, we add it to the final stimulus
toneparam.amplitude=amplitude;
toneparam.frequency=freqs(1);
toneparam=CalibrateTone(toneparam, cal);

tone=MakeTone(toneparam,samplerate);
sampledtoneduration=length(tone);
tones=zeros(nfreqs,sampledtoneduration);
tones(1,:)=tone;
for k=2:nfreqs
    toneparam.frequency=freqs(k);
    toneparam.amplitude=amplitude;
    toneparam=CalibrateTone(toneparam, cal);
    tones(k,:)=MakeTone(toneparam,samplerate);
end

%randomize phase %mw 11.19.2013
for k=1:nfreqs
    tones(k,:)=circshift(tones(k,:)', randi(ceil(samplerate/freqs(k))));
end

chord=sum(tones);
if ramp>0
    [edge,ledge]=MakeEdge(ramp,samplerate);         % and add the edge
    chord(1:ledge)=chord(1:ledge).*fliplr(edge);
    chord((end-ledge+1):end)=chord((end-ledge+1):end).*edge;
end

chord_starts=[0:(nchords-1)]';

sampled_start=max(1,sampled_start); % if sampled_start==0, we would have problems with indices below
chord_starts=sampled_start+chord_starts*(sampled_isi);

widths=0:sampledtoneduration-1;

idx=chord_starts(:,ones(1,sampledtoneduration))+widths(ones(1,nchords),:);
chord=chord(ones(1,nchords),:);

chordtrain(idx)=chord;

% !!!!! NOTE
% for now we're not checking the resulting amplitude!!! and we assume that
% the user is sane and able to check it himself/herself...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimparams=CalibrateTone(stimparams, cal);
findex=find(cal.logspacedfreqs<=stimparams.frequency, 1, 'last');
atten=cal.atten(findex);
stimparams.amplitude=stimparams.amplitude-atten;
