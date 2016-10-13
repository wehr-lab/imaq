global pref
Prefs('jls')    
cd('C:\Program Files\Teledyne DALSA\Sapera\CamFiles\User')
    camfilename='pantera_2.ccf';
    vid = videoinput('dalsa', 1, camfilename);

    src = getselectedsource(vid);
    %imaqmem(1e12);
    vid.timeout=5;
    vid.LoggingMode = 'memory';
    vid.FramesPerTrigger = 1;
vid.TriggerRepeat = Inf;
triggerconfig(vid,'hardware','risingEdge-ttl','trigger1');

vid.FramesPerTrigger = 10;
triggerconfig(vid,'immediate');

start(vid)
vid.FramesAcquired
stop(vid)

z = getframe(vid,9);