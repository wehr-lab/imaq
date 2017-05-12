function [filename,path]=MakeTuningCurve(numfreqs, minfreq, maxfreq, numamplitudes, ...
    minamplitude, maxamplitude, durations, ramp, include_whitenoise, isi, nrepeats)
%usage: [filename,path]=MakeTuningCurve(numfreqs, minfreq, maxfreq, ...
% numamplitudes, minamplitude, maxamplitude, duration, ramp, include_whitenoise, isi, nrepeats)
%creates an exper2 stimulus protocol file for a tuning curve stimulus
% inputs:
% numfreqs: number of frequency steps, log spaced between minfreq and maxfreq
%           use 0 for no tones (whitenoise only)
% minfreq: lowest frequency in Hz
% maxfreq: highest frequency in Hz
% numamplitudes: number of amplitude steps
% minamplitude: maximum amplitude in dB SPL (requires system to be calibrated)
% maxamplitude: maximum amplitude in dB SPL (requires system to be calibrated)
% durations: vector of different tone durations (in ms) (can be a single duration)
% ramp: on-off ramp duration in ms
% include_whitenoise: 0 or 1 to include white noise bursts at each amplitude
% isi: inter stimulus interval (onset-to-onset) in ms
% nrepeats: number of repetitions (different pseudorandom orders)
% outputs:
% creates a suitably named stimulus protocol in exper2.2\protocols
%
%
%example call: MakeTuningCurve(16, 1000, 32000, 3, 50, 80, 200, 10, 1, 500, 10)
%
%example call with multiple durations: 
%MakeTuningCurve(16, 1000, 32000, 3, 50, 80, [200 400],10,1, 500, 10) 
%
%this is good for a rate-level function (with 1 tone+WN):
%MakeTuningCurve(1, 7.3e3, 7.3e3, 6, 30, 80, [25], 5, 1, 500, 10)
%
%this is good for a rate-level function (with WN only):
%MakeTuningCurve(0, 0, 0, 9, 0, 80, 25, 3, 1, 500, 20)
% 
% This is what we've been using for binaural (with ear pieces in situ) mak 27Jul2010
% MakeTuningCurve(20, 1e3, 16e3, 5, 0, 70, [25], 3, 1, 350, 5)
% 
% Returns filename & path for chaining together sound, laser, vc makeprotocol functions.
% AKH 7/7/14

numdurations=length(durations);
logspacedfreqs = logspace( log10(minfreq) , log10(maxfreq) , numfreqs );
linspacedamplitudes = linspace( minamplitude , maxamplitude , numamplitudes );

if numfreqs==0; logspacedfreqs=[]; end

if include_whitenoise==1
    logspacedfreqs=[logspacedfreqs -1]; %add whitenoise as extra freq=-1
    numfreqs=numfreqs+1;
end

[Amplitudes,Freqs, Durations]=meshgrid( linspacedamplitudes , logspacedfreqs, durations );
neworder=randperm( numfreqs * numamplitudes * numdurations);
amplitudes=zeros(size(neworder*nrepeats));
freqs=zeros(size(neworder*nrepeats));
durs=zeros(size(neworder*nrepeats));

tdur=numfreqs * numamplitudes*numdurations *(mean(durations)+isi)/1000;%approx. duration per repeat

for nn=1:nrepeats
    neworder=randperm( numfreqs * numamplitudes * numdurations);
    amplitudes( prod(size(Amplitudes))*(nn-1) + (1:prod(size(Amplitudes))) ) = Amplitudes( neworder );
    freqs( prod(size(Freqs))*(nn-1) + (1:prod(size(Freqs))) ) = Freqs( neworder );
    durs( prod(size(Durations))*(nn-1) + (1:prod(size(Durations))) ) = Durations( neworder );
end

durstring=sprintf('%d-', durations);durstring=durstring(1:end-1);
%put into stimuli structure
stimuli(1).type='exper2 stimulus protocol';
if include_whitenoise
    stimuli(1).param.name= sprintf('Tuning curve +WN, %df(%d-%dHz)/%da(%d-%ddB)/%dd(%sms)/%dmsisi',...
        numfreqs,minfreq, maxfreq, numamplitudes,minamplitude, maxamplitude, numdurations, durstring,isi);
    stimuli(1).param.description=...
        sprintf('tuning curve, Tones +whitenoise, %d freq. (%d-%dkHz), %d ampl. (%d-%d dB SPL), %d durations (%sms), %dms ramp, %d repeats, %dms isi, %ds duration per repeat',...
        numfreqs, minfreq, maxfreq, numamplitudes,minamplitude, maxamplitude, numdurations, durstring, ramp, nrepeats, isi, round(tdur));
    filename=sprintf('tuning-curve-tones+WN-%df_%d-%dHz-%da_%d-%ddB-%dd_%sms-isi%dms-n%d',...
        numfreqs,minfreq, maxfreq, numamplitudes,minamplitude, maxamplitude, numdurations, durstring, isi, nrepeats);
else
    stimuli(1).param.name= sprintf('Tuning curve, %df(%d-%dHz)/%da(%d-%ddB)/%dd(%sms)/%dmsisi', ...
        numfreqs,minfreq, maxfreq, numamplitudes,minamplitude, maxamplitude, numdurations, durstring,isi);
    stimuli(1).param.description=sprintf('tuning curve, Tones only, %d freq. (%d-%dkHz), %d ampl. (%d-%d dB SPL), %d durations (%sms), %dms ramp, %d repeats, %dms isi, %ds duration per repeat',...
        numfreqs, minfreq, maxfreq, numamplitudes,minamplitude, maxamplitude, numdurations, durstring, ramp, nrepeats, isi, round(tdur));
    filename=sprintf('tuning-curve-tones-%df_%d-%dHz-%da_%d-%ddB-%dd_%sms-isi%dms-n%d',...
        numfreqs,minfreq, maxfreq, numamplitudes,minamplitude, maxamplitude, numdurations, durstring, isi, nrepeats);
end
for nn=1:length(amplitudes)
    if freqs(nn)==-1
        stimuli(nn+1).type='whitenoise'; %use nn+1 because stimuli(1) is name/description
        stimuli(nn+1).param.amplitude=amplitudes(nn);
        stimuli(nn+1).param.duration=durs(nn);
        stimuli(nn+1).param.ramp=ramp;
        stimuli(nn+1).param.next=isi;
    else
        stimuli(nn+1).type='tone';
        stimuli(nn+1).param.frequency=freqs(nn);
        stimuli(nn+1).param.amplitude=amplitudes(nn);
        stimuli(nn+1).param.duration=durs(nn);
        stimuli(nn+1).param.ramp=ramp;
        stimuli(nn+1).param.next=isi;
    end
end
global pref
Prefs
cd(pref.protocols)
save(filename, 'stimuli')
path=cd;


% keyboard