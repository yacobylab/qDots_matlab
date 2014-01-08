function gain=calibgain(CS, grpname, err) 

global fbdata; 
load('dBz_R_6626') 
fitfn=@(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(4))).*exp(-((x+p(6))/p(5)).^decay);
beta0=[0.015 0.008 10 0 50 0];
pars = fitwrap('',1:(size(data{1},2)),nanmean(data{1}(:),beta0,fitfn, [1 1 1 1 1 0]);
amp=pars(2); 


xval=plsinfo('xval',grpname,[],now); 
shottime=xval(2,1); 
if ~exist('amp','var')
    amp=0.01; 
end
freq=fbdata.gradtarget(2)./1e3; % dBz in GHz 
% p(1)=amp (V) p(2)=f (GHz) x=time (ns)
slpfn=@(p,x) 1e9*2*pi*p(1)*x*cos(2*pi*p(2)*x);  %this gives the slope in dv/df = V/Hz. 
slp=slpfn([amp freq],shottime);

Ns=CS.NSamples.get(); 
Sc=CS.SampleCount.get(); 

Gtrans=@(x) x.*2^(34)./(500*Ns.*Sc); 

sens=fbdata.VCOsens; 

gopt=1./(slp*sens); %The optimal gain
gain=Gtrans(gopt); 

dtheta=err/1e3*shottime; 
if abs(dtheta)>pi/2 
    warning('The error is greater than pi/2, Try a shorter time \n'
end

