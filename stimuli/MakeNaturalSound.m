function sound=MakeNaturalSound(param,fs)
% Return samples from soundfile referenced in params
% param should have
%   'file'      : pref.protocols/<path to .mat file>
%   'duration'  : duration of sound in ms
%   'amplitude' : the amplitude of the sound (0,1)
%   'next'      : ISI in ms
% fs = sampling rate in hz

% Should only be called once we have Pref'd
global pref

% Switch path seps if protocol made on different OS
if length(strsplit(param.file,filesep)) == 1
    if length(strsplit(param.file,'/')) > 1
        param.file = join(split(param.file,'/'),filesep);
    elseif length(strsplit(param.file,'\')) > 1
        param.file = join(split(param.file,'\'),filesep);
    end
end

cd(pref.protocols)

sound = load(char(param.file));
sound.sample.param.fs = 192000;
sample_fs = sound.sample.param.fs;
sound = sound.sample.sample;
sound=resample(sound, fs , sample_fs);

if length(sound)>(fs*param.duration/1000)
    sound = sound(1:(fs*param.duration/1000));
elseif length(sound)<(fs*param.duration/1000)
    sound(end+1:(fs*param.duration/1000)) = 0;
end


% For now we assume that amplitude has been adjusted during the protocol
% make stage. To change later -JLS101316




