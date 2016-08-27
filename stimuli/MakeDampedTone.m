function damped=MakeDampedTone(varargin)

global pref

% NOTE: this was the original header
% function damped=MakeDampedTone(Fc, halflife, duration, samplerate)

% [This function generates ramped and damped sinusoidal stimuli]
%
% function [ramped, damped]=rdamped(cf, halflife, period, stimulus_length, SR)
% Fc                -   carrier frequency (kHz)
% halflife          -   half-life (ms)
% period            -   period (ms)
% duration          -   duration (ms)
% samplerate        -   sampling rate (Hz)
% amplitude         -   sound pressure level of the sound (dB)

% orig: rdamped.m Wang Lab, Johns Hopkins University (Edited on January 13, 2004 by Tom Lu)

damped=[];

if nargin<2
    return;
end

params      =varargin{1};
samplerate  =varargin{2};
Fc          =params.frequency;
halflife    =params.halflife;
duration    =params.duration;
amplitude   =params.amplitude;

    amplitude=10*(10.^((amplitude-pref.maxSPL)/20));

    halflife=halflife/1000;      %halflife is now in sec

    npts=samplerate*(duration/1000);
    t=(1:npts)/samplerate;          % 'time' axis
    x = sin(2*pi*Fc*t)';            % carrier

    a = 1;
    b = 2;
    c = log(2);

    env = a*exp(-c*t/halflife)';    % envelope
    env = env/max(env);
%    env = env(1:length(x));
    env(end)=0;
    y=env .* x;
    y=amplitude.y;    
    damped = y;

% figure
% subplot(211)
% plot((1:length(damped))/samplerate*1000, damped)
% xlabel('Time (ms)')
% ylabel('amplitude')
% axis([0 length(damped)/samplerate*1000 -1 1])
% 
% subplot(212)
% plot((1:length(fft(damped)))/length(fft(damped))*(samplerate/1000), abs(fft(damped)))
% xlabel('Frequency (kHz)')
% ylabel('amplitude')
% axis([0 samplerate/1000 0 max(abs(fft(damped)))])