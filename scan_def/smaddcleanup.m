% function scan=smaddcleanup(scan, chan)
% Add (or edit) a cleanup function that restores (chan) to their
% value at the start of the scan.  

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

function scan=smaddcleanup(scan, chan)
  v=smget(chan);
  v=[v(:)];
  scan.cleanupfn(end+1).fn=@smconfigwrap;
  scan.cleanupfn(end).args={@smset,chan,v};
return
