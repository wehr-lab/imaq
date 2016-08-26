    cd('C:\Program Files\Teledyne DALSA\Sapera\CamFiles\User')
    camfilename='pantera.ccf';
    vid = videoinput('dalsa', 1, camfilename);

    src = getselectedsource(vid);
    imaqmem(1e12);
    vid.timeout=60;
    vid.LoggingMode = 'memory';
