function [snd,Fs] = pa_genripple3(vel,dens,md,durrip,durstat,varargin)
% [SND,FS] = PA_GENRIPPLE(VEL,DENS,MOD,DURDYN,DURSTAT)
%
% Generate a ripple stimulus with velocity (amplitude-modulation) VEL (Hz),
% density (frequency-modulation) DENS (cyc/oct), and a modulation depth MOD
% (0-1). Duration of the ripple stimulus is DURSTAT+DURRIP (ms), with the first
% DURSTAT ms no modulation occurring.
%
% These stimuli are parametrized naturalistic, speech-like sounds, whose
% envelope changes in time and frequency. Useful to determine
% spectro-temporal receptive fields.  Many scientists use speech as
% stimuli (in neuroimaging and psychofysical experiments), but as they are
% not parametrized, they are basically just using random stimulation (with
% random sentences).  Moving ripples are a complete set of orthonormal
% basis functions for the spectrogram.
%
% See also PA_GENGWN, PA_WRITEWAV

% 2011 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com
%
% Acknowledgements:
% Original script from Huib Versnel and Rob van der Willigen
% Original fft-script from NLS tools (Powen Ru)

%% Initialization
if nargin<1
	vel = 0; % (Hz)
end
if nargin<2
	dens = 1.2; % (cyc/oct)
end
if nargin<3
	md = 100; % Percentage (0-100%)
end
if nargin<4
	durrip = 1000; %msec
end
if nargin<5
	durstat = 500; %msec
end

%% Optional arguments
disp       = pa_keyval('display',varargin);
if isempty(disp)
	disp	= 1;
end
sv         = pa_keyval('save',varargin);
if isempty(sv)
	sv	= 'n';
end
plee       = pa_keyval('play',varargin);
if isempty(plee)
	plee	= 'n';
end
Fs         = pa_keyval('freq',varargin);
if isempty(Fs)
	Fs		= 48828.125; % Freq (Hz)
end
meth        = pa_keyval('method',varargin);
if isempty(meth)
	meth			= 'fft';
	% other meth = 'fft' works slightly faster
end
if strcmp(meth,'fft') % fft method is not defined for negative densities
	if dens<0
		dens	= -dens;
		vel		= -vel;
	end
end
md			= md/100; % Gain (0-1)

%% According to Depireux et al. (2001)
nFreq	= 128;
FreqNr	= 0:1:nFreq-1;
F0		= 250;
df		= 1/20;
Freq	= F0 * 2.^(FreqNr*df);
Oct		= FreqNr/20;                   % octaves above the ground frequency
Phi		= pi - 2*pi*rand(1,nFreq); % random phase
Phi(1)	= pi/2; % set first to 0.5*pi

%% Sounds
switch meth
	case 'time'
		rip		= genrip(vel,dens,md,durrip,Fs,nFreq,Freq,Oct,Phi);
	case 'fft'
		rip = genripfft(vel,dens,md,durrip,Fs,df);
end
if durstat>0 % if required, construct a static part of the noise, and prepend
	stat	= genstat(durstat,Fs,nFreq,Freq,Phi);
	% Normalize ripple power to static power
	nStat		= numel(stat);
	rms_stat	= norm(stat/sqrt(nStat));
	rms_rip		= norm(rip/sqrt(numel(rip)));
	ratio		= rms_stat/rms_rip;
	snd			= [stat ratio*rip];
else
	nStat = 0;
	snd = rip;
end



%% Normalization
% Because max amplitude should not exceed 1
% So set max amplitude ~= 0.8 (44/55)
snd			= snd/55; % in 3 sets of 500 trials, mx was 3 x 44+-1


%% Graphics
if disp
	plotspec(snd,Fs,nStat,Freq);
end

%% Play
if strcmpi(plee,'y');
	snd = pa_envelope(snd',round(10*Fs/1000));
	p	= audioplayer(snd,Fs);
	playblocking(p);
end

%% Save
if strcmpi(sv,'y');
	wavfname = ['V' num2str(vel) 'D' num2str(dens) '.wav'];
	pa_writewav(snd,wavfname);
end

function plotspec(snd,Fs,nStat,Freq)
close all;
t = (1:length(snd))/Fs*1000;
subplot(221)
plot(t,snd,'k-')
ylabel('Amplitude (au)');
ylabel('Time (ms)');
xlim([min(t) max(t)]);
hold on

subplot(224)
pa_getpower(snd(nStat+1:end),Fs,'orientation','x');
set(gca,'XTick',[0.5 1 2 4 8 16]*1000);
set(gca,'XTickLabel',[0.5 1 2 4 8 16]);
set(gca,'YTick',-160:20:0);
xlim([0.8*min(Freq) 1.2*max(Freq)])
ax = axis;
% xlim(0.6*ax([1 2]));
% set(gca,'YScale','log');
ylim([-160 0]);

subplot(223)
nsamples	= length(snd);
t			= nsamples/Fs*1000;
dt			= 12.5;
nseg		= t/dt;
segsamples	= round(nsamples/nseg); % 12.5 ms * 50 kHz = 625 samples
noverlap	= round(0.6*segsamples); % 1/3 overlap
window		= segsamples+noverlap; % window size
nfft		= 1000;
spectrogram(snd,window,noverlap,nfft,Fs,'yaxis');
% colorbar
cax = caxis;
caxis([0.7*cax(1) 1.1*cax(2)])
ylim([min(Freq) max(Freq)])
set(gca,'YTick',[0.5 1 2 4 8 16]*1000);
set(gca,'YTickLabel',[0.5 1 2 4 8 16]);
set(gca,'YScale','log');
drawnow

function snd = genstat(durstat,Fs,nFreq,Freq,Phi)
nTime	= round( (durstat/1000)*Fs ); % # Samples for Static Noise
time	= ((1:nTime)-1)/Fs; % Time (sec)


%% Modulate carrier, with static and dynamic part
snd		= 0;
for ii = 1:nFreq
	stat	= 1.*sin(2*pi* Freq(ii) .* time + Phi(ii));
	snd		= snd+stat;
end

function snd = genrip(vel,dens,md,durrip,Fs,nFreq,Freq,Oct,Phi)
nTime   = round( (durrip/1000)*Fs ); % # Samples for Rippled Noise
time	= ((1:nTime)-1)/Fs; % Time (sec)

%% Generating the ripple
% Create amplitude modulations completely dynamic without static
T			= 2*pi*vel*time;
F			= 2*pi*dens*Oct;
[T,F]		= meshgrid(T,F);
A			= 1+md*sin(T+F);
A			= A';

%% Modulate carrier
snd		= 0; % 'initialization'
for ii = 1:nFreq
	rip		= A(:,ii)'.*sin(2*pi*Freq(ii) .* time + Phi(ii));
	snd		= snd+rip;
end

function snd = genripfft(vel,dens,md,durrip,Fs,df)
% Generate ripple in frequency domain

%%
Ph		= 0-pi/2;

%% excitation condition
durrip	= durrip/1000;	% ripple duration (s)
f0		= 250;	% lowest freq
BW		= 6.4;	% bandwidth, # of octaves
Fs		= round(Fs*2)/2; % should be even number
maxRVel	= max(abs(vel(:)),1/durrip);

%% Time axis
tStep		= 1/(4*maxRVel);
tEnvSize	= round((durrip/tStep)/2)*2; % guarantee even number for fft components
tEnv		= (0:tEnvSize-1)*tStep;

%% Frequency axis
oct			= (0:round(BW/df*2)/2-1)*df;
fr			= pa_oct2bw(f0,oct)';
fEnv		= log2(fr./f0);
fEnvSize	= length(fr);	% # of component

%% Compute the envelope profile
ripPhase	= Ph+pi/2;
fPhase		= 2*pi*dens*fEnv + ripPhase;
tPhase		= 2*pi*vel*tEnv;
A			= md*(sin(fPhase)*cos(tPhase)+cos(fPhase)*sin(tPhase));
A			= 1+A; % shift so background = 1 & profile is envelope

%% freq-domain AM
nTime		= durrip*Fs; % signal time (samples)

%% roll-off and phase relation
th			= 2*pi*rand(fEnvSize,1); % component phase, theta
S			= zeros(1, nTime); % memory allocation
tEnvSize2	= tEnvSize/2;
for ii = 1:fEnvSize
	f_ind				= round(fr(ii)*durrip);
	S_tmpA				= fftshift(fft(A(ii,:)))*exp(1i*th(ii))/tEnvSize*nTime/2;
	
	pad0left		= f_ind - tEnvSize2 - 1;
	pad0right		= nTime/2 - f_ind - tEnvSize2;
	if ((pad0left > 0) && (pad0right > 0) )
		S_tmpB = [zeros(1,pad0left),S_tmpA,zeros(1,pad0right)];
	elseif ((pad0left <= 0) && (pad0right > 0) )
		S_tmpB = [S_tmpA(1 - pad0left:end),zeros(1,pad0right)];
	elseif ((pad0left > 0) && (pad0right <= 0) )
		S_tmpB = [zeros(1,pad0left),S_tmpA(1:end+pad0right)];
	end
	S_tmpC	= [0, S_tmpB, 0, fliplr(conj(S_tmpB))];
	S		= S + S_tmpC; % don't really have to do it all--know from padzeros which ones to do...
end
snd = real(ifft(S));

