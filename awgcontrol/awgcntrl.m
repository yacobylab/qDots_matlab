function awgcntrl(cntrl, chans)
% awgcntrl(cntrl, chans)
% cntrl: stop, start, on off, wait, raw|amp, extoff|exton

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global awgdata;

if nargin <2
    chans = awgdata.chans;
end

switch cntrl
    case 'stop'
        fprintf(awgdata.awg, 'AWGC:STOP');
        
    case 'start'
        fprintf(awgdata.awg, 'AWGC:RUN');
        
    case 'off'
        for i = chans
            fprintf(awgdata.awg, 'OUTPUT%i:STAT 0', i);
        end
        
    case 'on'
        for i = chans
            fprintf(awgdata.awg, 'OUTPUT%i:STAT 1', i);
        end
        
    case 'wait'
        to = awgdata.awg.timeout;
        awgdata.awg.timeout = 600;
        query(awgdata.awg, '*OPC?');
        awgdata.awg.timeout = to;
        
    case 'raw'
        %awgcntrl('stop');
        %awgcntrl('off');
        for i = chans
            fprintf(awgdata.awg, 'AWGC:DOUT%i:STAT 1', i);
        end
        %awgcntrl('on');
        %awgcntrl('start');
        
    case 'amp'
        %awgcntrl('stop');
        %awgcntrl('off');
        for i = chans
            fprintf(awgdata.awg, 'AWGC:DOUT%i:STAT 0', i);
        end
        %awgcntrl('on');
        %awgcntrl('start');
        
    case 'exton'    %adds external DC to outputs specified in chans
        %awgcntrl('stop');
        %awgcntrl('off');
        for i = chans
            fprintf(awgdata.awg, 'SOUR%i:COMB:FEED "ESIG"', i);
        end
        %awgcntrl('on');
        %awgcntrl('start');
        
   
    case 'extoff'   %turns off external DC
        %awgcntrl('stop');
        %awgcntrl('off');
        for i = chans
            fprintf(awgdata.awg, 'SOUR%i:COMB:FEED ""', i);
        end
        %awgcntrl('on');
        %awgcntrl('start');
        
        
end
