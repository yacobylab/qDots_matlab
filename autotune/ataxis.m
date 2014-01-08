% Swap the X and Y axes on tunedata.  Need to manually change the pulse
% group

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

function ataxis(xchan, ychan, rchan)
global tunedata;
  if(~ischar(xchan))
      xchan=sprintf('%d',xchan);
  end
  if(~ischar(ychan))
      ychan=sprintf('%d',ychan);
  end
  
  tunedata.chrg.scan.cleanupfn(1).args{2}{1}(end)=xchan;
  tunedata.chrg.scan.cleanupfn(1).args{2}{4}(end)=xchan;
  tunedata.chrg.scan.cleanupfn(1).args{2}{2}(end)=ychan;
  tunedata.chrg.scan.cleanupfn(1).args{2}{3}(end)=ychan;
  tunedata.chrg.scan.loops(1).setchan{1}(end) = xchan;
  tunedata.chrg.scan.loops(2).setchan{1}(end) = ychan;
  tunedata.lead.scan.cleanupfn = tunedata.chrg.scan.cleanupfn;
  tunedata.lead.scan.loops(2).setchan{1}(end) = xchan;
  tunedata.lead.scan.loops(2).setchan{2}(end) = ychan;
  tunedata.lead.scan.loops(2).getchan={rchan};
  
  tunedata.line.scan.cleanupfn = tunedata.chrg.scan.cleanupfn;
  tunedata.line.scan.loops(1).setchan{1}(end)=xchan;
  tunedata.line.scan.loops(1).setchan{2}(end)=ychan;
  tunedata.line.scan.loops(2).getchan={rchan};  
  
  tunedata.load.scan.cleanupfn = tunedata.chrg.scan.cleanupfn;
  tunedata.load.scan.consts(1).setchan(end)=xchan;
  tunedata.load.scan.consts(2).setchan(end)=ychan;
  tunedata.load.scan.loops(1).getchan{1}=rchan;
  
  tunedata.read.scan.cleanupfn = tunedata.chrg.scan.cleanupfn;
  tunedata.read.scan.loops(2).setchan{1}(end)=xchan;
  tunedata.read.scan.loops(2).setchan{2}(end)=ychan;  
  tunedata.read.scan.loops(2).getchan{1}=rchan;
  
  tunedata.tl.scan.consts(1).setchan(end)=xchan;
  tunedata.tl.scan.consts(2).setchan(end)=ychan;
  tunedata.tl.scan.loops(1).getchan{1}=rchan;
  
  tunedata.stp.scan.consts(1).setchan(end)=xchan;
  tunedata.stp.scan.consts(2).setchan(end)=ychan;
  tunedata.stp.scan.loops(1).getchans{1}=rchan;
end
