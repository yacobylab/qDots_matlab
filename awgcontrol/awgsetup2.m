
% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

global awgdata;

if exist('smdata', 'var') && isfield(smdata, 'inst')
    awgdata.awg = smdata.inst(sminstlookup('AWG5000')).data.inst;
end
if strcmp(computer, 'GLNX86')
    %awgdata.datafile = '~ygroup/qDots/pulsedata_10mV_0708.mat';
    %awgdata.datadir = '~ygroup/qDots/awg_pulses/nfs/';
    awgdata.datafile = '~ygroup/qDots/awg_pulses/pulsedata_1111.mat';
    awgdata.grpdir = '~ygroup/qDots/awg_pulses/sequence_files/';
    awgdata.datadir = '~ygroup/qDots/awg_pulses/waveforms/2010_11_09/';
else
    awgdata.datafile = 'z:/qDots/awg_pulses/pulsedata_1111.mat';
    awgdata.grpdir = 'z:/qDots/awg_pulses/sequence_files/';
    awgdata.datadir = 'y:/2010_11_09/';    
end
awgdata.chans = 1:4; % default channels for control operations
awgdata.nchan = 4;
awgdata.offp = 1; % ref pulse in order to include a blank. Seems to help with noise.
awgdata.trigp = 2;
awgdata.plsgen = 'A';
awgreadxval;
awgsyncdata('load');
awgcalcshift(1:length(awgdata.pulsedata));
%%
if 0

    smprintf(10, 'MMEM:MSIS "C:"')
    smprintf(10, 'MMEM:CDIR "\Documents and Settings\OEM\My Documents\waveforms\2010_11_09"')
    % 10 +33 dB atten, 10 mV
    for i = 1:4
        smprintf(10, 'SOUR%i:VOLT:AMPL .6', i);
    end
    %smprintf(16, 'SOUR1:VOLT:AMPL 1.405');
    %smprintf(16, 'SOUR2:VOLT:AMPL 1.385');
    
end
