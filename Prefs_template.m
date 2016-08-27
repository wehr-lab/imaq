function Prefs
% sets the standard exper preferences
% should be called at the very beginning of the ExperSomeName.m:-)


% 'global pref' should be added to the beginning of each module
global pref

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT USER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pref.username contains the default username; it's actually used only by
% modified Control module to stop bothering the user with that damn empty
% dialog box in the beginning
% :( the dialog box is back again mw 051209
options.WindowStyle='modal';
if ~isfield(pref, 'username')
    pref.username=inputdlg('Please enter your username','username',1,{'lab'}, options);
    if ~isempty(pref.username); pref.username=pref.username{:};end
end
while isempty(pref.username)
    pref.username=inputdlg('Please enter your username','username',1,{'lab'}, options);
    if ~isempty(pref.username); pref.username=pref.username{:};end
end
pref.rig='rig2';
pref.usebak=0;
pref.mkdir=0;

pref.numchannels=2; %flag for whether to process a second data channel in E2ProcessDAQFile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pref.home stores path to the user's home directory
pref.home='c:\lab\';

% pref.processed_data_dir stores path to the user's processed_data directory
pref.processed_data_dir=fullfile(pref.home, sprintf('Data-%s-processed', pref.username));

% pref.experhome stores path to the exper's main tree
pref.experhome=[pref.home 'exper2.2\'];

% paths to three main directories inside the exper tree
% pref.modules is the directory with all available exper modules
pref.modules=[pref.experhome 'modules\'];
% pref.utility contains helper functions
pref.utility=[pref.experhome 'utility\'];
% pref.rp2 contains 'exper RP2 circuit file' (see pref.rp2soundcircuit below) and 
% utilities for generating sound stimuli 
pref.tdt    =[pref.experhome 'rp2\'];
% pref.stimulifunctions contains the path to the dir with the stimuli
% generating stuff
pref.stimulifunctions=[pref.experhome 'stimuli\'];

%pref.patchpreprocessamp sets default for patchpreprocess "amp" param
%that sets whether the gain is read from the axopatch or set manually
if strcmp(pref.username,'ira')
    pref.patchpreprocessamp=1; %1 for "use axopatch", 2 for "manual gain"
else
    pref.patchpreprocessamp=2; %1 for "use axopatch", 2 for "manual gain"
end

% path to default data directory
pref.data=[pref.home sprintf('data-%s\\', pref.username)];
% path to default stimuli directory
pref.stimuli=[pref.experhome 'protocols\'];

% path to default directory with behavior protocols for BoxMaster
pref.boxprotocols=[pref.experhome 'boxprotocols\'];

%Drive Letter for Open Ephys Data 
pref.openephysdrive='d';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HARDWARE CHANNELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% DEFAULT MODULES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pref.default_modules contains a list of modules that will be open when
% you start Exper2
needed={'ai';...               %analog input
        'ao';...               %analog output
        'dio';...              %digital input/output
        'patchpreprocess';...  %axopatch settings and channel info
        'stimulusprotocol';...
        'pathdisplay';...
        };
% The most current module list is in the control module after you runexper
% >> control >> Modules 
if strcmp(pref.username,'mak')
    wanted={'calibratespeaker';...
            'ppalaser';...
            };
    % elseif strcmp(pref.username,'add your name here')
else % standard modules for most experiments
    wanted={'sealtest';...
        'onliner';...
        'sealtest';...
        'holdcmd';...
        };
end

pref.default_modules=[needed;wanted];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RP2 STUFF:-)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default setings for the TDT system
% For each TDT thing (RP2, RM1,...) pref.rp2 contains:
% type: usualy RP2, we might also use RM1 (MUST BE IN CAPITAL LETTERS!!!)
% interface: how we are connected to the TDT stuff; either GB (Gigabit) or
% USB (guess what)
% circuit: filename of the circuit that will be loaded to the RP2/RM1/...
% description: circuit function. Modules using RP2 can find out which RP2
% to use
% IMPORTANT NOTES: 
% The first RP2/RM1 is the main one by default, which means that it will be
% used to produce the stimuli, i.e. it'd better be connected to the speaker
% and you'd better specify the correct circuit.
% The circuit files MUST be placed in pref.tdt directory (see above)
%      Type  Interface      Circuit                   Description
rp2={  'RP2'   'USB'    'soundrp2-1speaker.rco'      'sound';...
    };
% rp2={  'RP2'   'USB'    'soundrp4-simple-self-falling-edge-final.rco'                   'sound';...
%        'RP2'   'USB'    'timestamp-22bit-and-water-delivery-falling-edge-self.rco'      'behavior, timestamp';...
%     };

for n=1:size(rp2,1)
    pref.rp2(n).type         =rp2{n,1};
    pref.rp2(n).interface    =rp2{n,2};
    pref.rp2(n).circuit      =rp2{n,3};
    pref.rp2(n).description  =rp2{n,4};
end

% name of the RP2 sound circuit file
% the file must be placed in the pref.tdt directory (see above)
pref.ppfilterfile=[pref.tdt 'PPfilter.mat'];

% max. sound amplitude
% RP2 can deliver sound with amplitudes in range of +/-10V. The max.
% amplitude corresponds to some sound pressure level, which depends on the
% speaker, equalizer, etc.
%pref.maxSPL=71;
pref.maxSPL=100;

% RP2 sampling rate we use
pref.RP2Fs=50e6/512;

% description of the various speakers connected to rp2s
%           RP2#    Channel#
speakers={   1          1;...
             1          2;...
         };
     
for n=1:size(speakers,1);
    pref.speakers(n).rp2    =speakers{n,1};
    pref.speakers(n).channel=speakers{n,2};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOUNDCARD STUFF :-0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of soundcard output channels you want to use 
%cannot exceed capability of the soundcard you have installed
%pref.num_soundcard_outputchannels=3;
pref.num_soundcard_outputchannels=4; %mw 081412 for PPA Laser channel 
%help PPALaser for more info

%requested latency class for PPA initialization 
%this differs by machine and the wrong choice can cause
%dropouts, for now we find the best value empirically  
pref.reqlatencyclass=1; %4 is best for rig1, but rig 2 seems to do better with 1, right??? 
     %We've now tried:0, 1, 10, (None work, we keep getting a double
     %trigger that leads to the error message when trying to PlotTC:
     %"??? Subscripted assignment dimension mismatch."
     %MAK,8Jul2009,9:52AM
%use printdevices.m to figure out which device id to use for your soundcard
%use the one called "ASIO Lynx"

%for PPASound
pref.SoundFs=192e3;
pref.soundmethod='PPAsound'; %choose from 'AOSound', 'PPAsound', or 'soundmachine'

pref.soundcarddeviceID=GetAsioLynxDevice;
%added new function GetAsioLynxDevice so we don't have to keep switching
%this setting around every time windows changes the ASIO Lynx deviceID.
%pref.soundcarddeviceID=26; %for PPASound; for some reason on the new rig 2 
%switched from 24 to 26 mw 01-24-2012
% cpu this number alternates between 22 & 24 mak15feb2011
% now it is 22 again mak 18feb2011. Perhaps Whitney's or Sam's experiments
% cause it to change???
% Use PrintDevices.m to find the number being currently used
% 05/28/11 printdevice ID switched to 22. changed back to 24 manually XG
% 05/30/11 printdevice ID switched to 24. changed back to 22 manually XG
% 10/18/11 printdevice ID switched to 24. changed back to 22 manually XG
%010512 changed from 26 to 24
%for AOSound
%pref.soundcarddeviceID=0; %for AOSound
%pref.soundmethod='AOSound'; %choose from 'AOSound',%soundmethod chooses which soundcard interface to use
%pref.SoundFs=96e3;
%don't forget to change the soundcarddeviceID too!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOW TO CREATE THE STIMULI
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEFORE WE FINISH...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the 'standard' paths must be put into the Matlab's path
addpath(pref.home);
addpath(pref.modules);
addpath(pref.utility);
addpath(pref.tdt);
addpath(pref.stimulifunctions);
addpath(pref.boxprotocols);

% set(0,'DefaultFigureColor',[.925 .914 .847]);
% set(0,'DefaultUicontrolBackgroundColor',[.925 .914 .847]);
% set(0,'DefaultuimenubackgroundColor',[.925 .914 .847]);


