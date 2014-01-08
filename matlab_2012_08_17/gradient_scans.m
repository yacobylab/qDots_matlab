%% Find gradients
% Note awesome code that auto-restores most recent channel values.
function [scan] = gradient_scans();
if 1
  sweepchans={'PlsRamp2','PlsRamp1'};
  stepchans={'1a','2a','3a','4a','1b','2b','3b','4b','N12','N34','T12','T34','PlsRamp3','PlsRamp4','PlsRamp2','PlsRamp1'};
  step=     [2   ,2   ,5   ,10   ,2   ,2   ,5   ,10   ,1    ,7    ,1    ,7    ,8       ,8       ,1       ,1   ]*1e-3;
else
  sweepchans={'PlsRamp3','PlsRamp4'};
  stepchans={'1a','2a','3a','4a','1b','2b','3b','4b','N12','N34','T12','T34','PlsRamp1','PlsRamp2','PlsRamp3','PlsRamp4'};
  step=     [10  ,5   ,2   ,2    ,10  ,5   ,2   ,2    ,7    ,1    ,7    ,1    ,8       ,8       ,1       ,1   ]*1e-3;    
%stepchans={'count','count2'};
end
rng=[-5e-3 5e-3];
gradbase = smget(stepchans);
gradbase = [gradbase{:}];
for i=1:length(sweepchans)
  scan(i).loops(1).setchan=sweepchans{i};
  scan(i).loops(1).rng = rng;
  scan(i).loops(1).npoints = 200;
  scan(i).loops(1).trigfn.fn=@smatrigAWG;
  scan(i).loops(1).trigfn.args={[16]};
  scan(i).loops(1).ramptime=-0.005;
  scan(i).loops(2).setchan=stepchans;
  for m=1:length(stepchans)
    scan(i).loops(2).trafofn(m).fn=@stepper;
    scan(i).loops(2).trafofn(m).args={m,step(m),gradbase(m)};    
  end
  scan(i).loops(2).rng = 0:(length(stepchans)*2-1);
  scan(i).loops(2).npoints = length(stepchans)*2;
  scan(i).loops(2).getchan={'DAQ1'};
  scan(i).loops(2).ramptime=[];
  scan(i).loops(3).npoints=3;
  scan(i).loops(3).setchan={'count'};
  scan(i).configfn(1).fn=@pre;
  scan(i).configfn(1).args={stepchans};
  scan(i).cleanupfn(1).fn=@cleanup
  scan(i).cleanupfn(1).args={stepchans,gradbase};
  scan(i).configfn(2).fn=@smabufconfig2;
  scan(i).configfn(2).args={'arm', [1]};
  scan(i).consts(1).setchan='samprate';
  scan(i).consts(1).val=8e6;
  scan(i).consts(2).setchan=sweepchans{mod(i,2)+1};
  scan(i).consts(2).val=scan(i).loops(1).rng(1);
  scan(i).disp(1).loop=2;
  scan(i).disp(1).channel=1;
  scan(i).disp(1).dim=1;
end

return;
function [scan]=pre(scan,stepchans)
  gradbase = smget(stepchans); 
  gradbase=[gradbase{:}];
  scan.cleanupfn(1).args{2}=gradbase;
  for l=1:length(gradbase)
    scan.loops(2).trafofn(l).args{3}=gradbase;
  end
return;
function [scan]=cleanup(scan,stepchans,values)
  smset(stepchans,values);
return
function val=stepper(x,y,i,s,gradbase)
  ch = floor(x(2)/2)+1;
  sgn = mod(x(2),2)-0.5;
  gradbase(ch)=gradbase(ch)+sgn*s;
  val=gradbase(i);
return;