function awginit

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global awgdata;

awgdata.awg = tcpip('yacobyawg.physics.harvard.edu', 4000);
awgdata.awg.ByteOrder = 'littleEndian';
awgdata.awg.timeout = 1;
awgdata.awg.OutputBufferSize = 2^18;
awgdata.datafile = 'z:/qDots/pulsedata_0308.mat';
