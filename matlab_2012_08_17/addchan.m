function scan=addchan(scan,chan)
%add a getchan chan (str) to the first loop of scan. fConfSeq2 puts the histogramming in the
%3rd data cell and does not allow more cells. This fcn changes the procfn to put the
%histogramming in a 4th cell, and adds chan as the 3rd getchan. 
scan.loops.procfn(4)=scan.loops.procfn(3); 
scan.loops.procfn(3).dim=[]; scan.loops.procfn(3).fn=[];
scan.loops.procfn(1).fn(1).outchan=4; 
scan.disp(2).loop=1; scan.disp(2).channel=3; scan.disp(2).dim=1; 
scan.disp(1).loop=1; scan.disp(1).channel=1; scan.disp(1).dim=2; 
scan.loops.getchan{3}=chan;