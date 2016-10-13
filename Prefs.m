function Prefs(user)
% sets the standard exper preferences
% should be called at the very beginning of the ExperSomeName.m:-)
% 'global pref' should be added to the beginning of each module
global pref

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIG SPECIFIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pref.ccf = 'C:\Program Files\Teledyne DALSA\Sapera\CamFiles\User\pantera_2.ccf';
pref.trigger_type = {'hardware','risingEdge-ttl','trigger1'}; % from triggerinfo on the vid object
pref.fps = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT USER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isa(user,'cell')
    user = user{1};
end
pref.username = user;
pref.rig      = 'rig2';
pref.usebak   = 0;
pref.mkdir    = 0;
% Move to 'hw settings' or something
pref.numchannels = 2; %flag for whether to process a second data channel in E2ProcessDAQFile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off MATLAB:MKDIR:DirectoryExists

pref.base='c:\lab\';
pref.home='c:\lab\imaq\';

% Anything that we can't guarantee exists we check and we make
% Declare dirs to make like 'name','subdir below pref.base'
alldirs = {
    'data',sprintf('data\\%s\\', pref.username);
    'processed_data',sprintf('data\\%s-processed\\', pref.username);
    'stimuli','stimuli\';
    'protocols','protocols\';
    'calibration','calibration\'
};

% Make alldirs and assign to prefs fields
for i = 1:size(alldirs,1)
   cd(pref.base)
   
   % If we were given nested folders, get to the right level
   dirsplit = strsplit(alldirs{i,2},'\');
   if length(dirsplit)>2
       for j=1:(length(dirsplit)-2); % We've been including trailing backslashes, which are returned as empty strings by strsplit
           if ~exist([pwd,dirsplit{j}],'dir')
               mkdir(dirsplit{j});
               cd(dirsplit{j});
           else
               cd(dirsplit{j});
           end
       end
   end
   
   if ~exist([pwd, dirsplit{end-1}],'dir')
       mkdir(dirsplit{end-1});
   end
   
   pref.(alldirs{i,1}) = [pref.base, alldirs{i,2}];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HARDWARE PREFS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pref.fs = 192000;
pref.dev_id = GetAsioLynxDevice;
pref.n_chan = 4;
pref.buff_size = 512;
pref.soundmethod='PPAsound'; %choose from 'AOSound', 'PPAsound', or 'soundmachine'
pref.maxSPL = 60;
pref.runMode = 0; % Turns off soundcard after playback





% Here you can specify AI channels (where the stuff gets in) and AO
% channels (where the stuff gets out)
% The format is: 
% Each hardware channel is described in one row of ai_channels or ao_channels 
% Each row contains the following:
% Channel# (integer): number of the channel in the data acquisition board (break-out box) - where you put the cable;
% Channel 'type' (string): characteristic string associated with channel - something that describes the channel 
% Channel name (string): human readable characteristics - put your favourite name here
% Channel status (string): either permanent or temporary. Permanent ones are used
% all the time, temporary only sometimes (to get some values from the Axopatch, for example)
% Channel color: only used for plotting data from various channels
%          channel#     channel 'type'          channel name           channel status      chanel color  
ai_channels={   0       'datachannel-patch'     'AxopatchData1'        'permanent'         [0 .5 0];...
                1       'stimulichannel'        'TDT'                  'permanent'         'k';...
                2       'triggerchannel'        'Triggers'             'permanent'         'b';...
                3       'datachannel2-patch'    'AxopatchData2'        'permanent'         'm';...
                4       'soundcardtrigchannel'  'soundcardtrigchannel' 'permanent'         'b';...
                5       'gainchannel'           'AxopatchGain1'        'temporary'         'b';...
                7       'modechannel'           'AxopatchMode1'        'temporary'         'b';...
                   };
               %4       'microphonechannel' 'Microphone'        'permanent';... 
        
 ao_channels={   0       'commandchannel'    'AxopatchCommand1'   'permanent'         'b';...
                1       'ledchannel'        'LEDChannel'         'permanent'         'g';...
%  ao_channels={   0       'ledchannel'        'LEDChannel'         'permanent'         'g';...
%                  1       'commandchannel'    'AxopatchCommand1'   'permanent'         'b';...
                 };
%                 1       'ttlchannel'        'TTLChannel'        'permanent'         'b';...

% in exper (the software, NOT the structure:-) you can find the AI (AO) channels in pref.ai_channels (pref.ao_channels)                              
for n=1:size(ai_channels,1)
    pref.ai_channels(n).number =ai_channels{n,1};
    pref.ai_channels(n).channel=ai_channels{n,2};
    pref.ai_channels(n).name   =ai_channels{n,3};
    pref.ai_channels(n).status =ai_channels{n,4};    
    pref.ai_channels(n).color  =ai_channels{n,5};    
end

for n=1:size(ao_channels,1)
    pref.ao_channels(n).number =ao_channels{n,1};
    pref.ao_channels(n).channel=ao_channels{n,2};
    pref.ao_channels(n).name   =ao_channels{n,3};
    pref.ao_channels(n).status =ao_channels{n,4};    
    pref.ao_channels(n).color  =ai_channels{n,5};    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIMULUS FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here you can set different stimulus types and names of the functions that
% create them. In general, all these functions expect any number of
% parameters. In reality, functions making various sounds should receive
% one parameter called params (see each function for details) and one
% called samplerate. For the other functions, see them. And, obviously,
% these functions should live somewhere in the Path. Usually they live in
% pref.stimulifunctions (see above).
%                      stimulus type    stimulus function           stimulus trigger type
pref.stimulitypes={     'tone'          'MakeTone'                  'sound';...
                        'bintone'       'MakeBinTone'               'sound';...
                    	'2tone'         'Make2Tone'                 'sound';...
                        'whitenoise'    'MakeWhiteNoise'            'sound';...
                        'binwhitenoise' 'MakeBinWhiteNoise'         'sound';...
                        'amnoise'       'MakeAMNoise'               'sound';...
                        'amtone'        'MakeAMTone'                'sound';...
                        'fmsweep'       'MakeFMSweep'               'sound';...
                        'fmtone'        'MakeFMTone'                'sound';...
                        'dampedtone'    'MakeDampedTone'            'sound';...
                        'rampedtone'    'MakeRampedTone'            'sound';...
                        'noise'         'MakeNoise'                 'sound';...
                        'spectrogram'   'MakeSoundFromSpectrogram'  'sound';...
                        'naturalsound'  ''                          'sound';...
                        'clicktrain'    'MakeClickTrain'            'sound';...
                        'chordtrain'    'MakeChordTrain'            'sound';...
                        'tonetrain'     'MakeToneTrain'             'sound';...
                        'oddball'       'MakeOddball'               'sound';...
                        'randomchords'  'MakeRandomChords'          'sound';...
                        'soundfile'     ''                          'sound';...
                        'ASR'           'MakeASR'                   'sound';...
                        'GPIAS'         'MakeGPIAS'                 'sound';...
                        'itdwhitenoise' 'MakeITDWhiteNoise'         'sound';...
                        'itdtone'       'MakeITDTone'               'sound';...
                        'gapinnoise'    'MakeGapInNoise'            'sound';...
                        % now something for pulses
                        'pulse'         'MakePulse'                 'ao';...
                        'aopulse'       'MakeAOPulse'               'ao';...
                        'holdcmd'       'MakeHoldCmd'               'ao';...
                        'led'           ''                          'ao';...
                        %LED output
                        
                        
                        % something for water delivery. We don't need any
                        % function, but it has to be here, because it is a
                        % valid stimulus type!
                        'waterdelivery' ''                          'waterdelivery';...
                  
                        % something for visual stimuli (psych toolbox)
                        % mw072108
                        'grating'       ''                          'visual'
                        };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OTHER STUFF...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here we can set the general pattern for expid
% each recording will add -xxx to expid to distiguish between different
% recordings
% For example I use: YYYYMMDD-'animalid'-'recording#' (ie
% 20030101-thr999-xxx), where -xxx is added by exper
% Animal id can be changed by ExperLog module, for example
expid=strrep(datestr(date,2), '/', '');
animalid=pref.username;

pref.expid=[expid '-' animalid];




