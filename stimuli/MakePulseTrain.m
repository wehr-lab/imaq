function clicktrain=MakePulseTrain(varargin)

global exper pref

% Creates a click train, but unlike MakeClickTrain which creates trains of
% short ramped bursts of white noise, this produces trains of short
% rectangular pulses. Noise bursts are fine for durations > 5ms or so, but
% for very short durations the finite sampling causes the actual bursts to
% vary quite a bit in spectrum and level. The pulses allow better control
% of duration and level for clicks shorter than ~1 ms.

% Input
%   nclicks         -   number of clicks
%   isi             -   inter-stimulus interval, i.e. interval between the
%                       start of previous click and start of the next click
%   clickduration   -   duration of an individual click (ms)
%   amplitude       -   click amplitude (dB)
%   start           -   start of the first click after the trigger (ms)
%   duration        -   total duration of the click train (ms) %%???not used? mw 04.13.06
%   next            -   inter-click-train-interval, i.e. when the next
%                       click train should follow the previous one (ms)
%   samplerate      -   sampling rate (Hz)
% Output
%   pulsetrain   -   the required square waveform
%
%modified from MakeClickTrain, mw 060811
    
clicktrain=[];

if nargin<2
    return;
end

params=varargin{1};
samplerate=varargin{2};
start           =params.start/1000;              % in s
clickduration   =params.clickduration/1000;      % in s
amplitude       =params.amplitude;
nclicks         =params.nclicks;
isi             =params.isi/1000;                % in s

    train_length=(start+clickduration+(nclicks-1)*isi); % in s
    
    %prepare the samples
    clicktrain=zeros(ceil(train_length*samplerate),1);
    
    sampled_duration=round(clickduration*samplerate);
    sampled_start=floor(start*samplerate);
    sampled_isi=round(isi*samplerate);

    click=ones(1,sampled_duration);       % corresponds to t=0:1/samplerate:duration;
% % %     [edge,ledge]=MakeEdge(ramp,samplerate);         % and add the edge
% % %     click(1:ledge)=click(1:ledge).*fliplr(edge);
% % %     click((end-ledge+1):end)=click((end-ledge+1):end).*edge;
    
    click_starts=[0:(nclicks-1)]';
    sampled_start=max(1,sampled_start); % if sampled_start==0, we would have problems with indices below
    click_starts=sampled_start+click_starts*(sampled_isi);
    
    widths=0:sampled_duration-1;
    
    idx=click_starts(:,ones(1,sampled_duration))+widths(ones(1,nclicks),:);
    click=click(ones(1,nclicks),:);
    
    clicktrain(idx)=click;
    
    clicktrain=clicktrain./(max(abs(clicktrain)));                 % normalize, so we could fit to +/-10V

%     amplitude=10*(10.^((amplitude-pref.maxSPL)/20));
    amplitude=1*(10.^((amplitude-pref.maxSPL)/20)); %mw 102908
    clicktrain=amplitude.*clicktrain;
