a= [];
for j = 1:10
    a(end+1) = sm_getgradient;
end

%%
epsvals=  [2.4:.2:4];

for j = epsvals
clear pg 
pg.pulses=55;
pg.ctrl='loop pack'; 
t1grp = {};
%Parameters: [pulselength, eps, time(us)] 

eps = j;
func = @plsupdate;

pg.dict={struct('prep','@null'),'right'};
pg.chan = [3 4]; 


pg.varpar=(linspace(1,100,64))';%128+0*(1:128)';
pg.params = [110, eps, 0];
pg.xval = eps;
pg.name=('T1_eps_gnd_2013_04_06_R');
t1grp=[t1grp, pg.name];
func(pg);

pg.name=('T1_eps_ex_2013_04_06_R');
pg.dict={struct('prep','@dbzpi'),'right'};
t1grp=[t1grp, pg.name];
func(pg); 
awgadd(t1grp);
awgcntrl('on start wait err raw')

if sm_setgradient
    fprintf('epsilon = %.2d\n',j);
    d=smrun(fConfSeq2([28:30],struct('nloop',5,'nrep',1000, 'opts','pol','datachan','DAQ2')),smnext('T1_eps_R'));
    if any(isnan(d{1}))
       break; 
    end
    if abs(sm_getgradient-fbdata.gradtgt(2))> 10
       fprintf('gradient likely came unlocked for eps = %.2d \n',j); 
    end
else
    fprintf('could not lock gradient for eps = %.2d\n',j);
end
end
sleep

%% echo
awgrm(27,'after'); awgclear('unused');
awgadd('dBz_12_64_1_R');
for vrf = [3.1, 3.3 3.5]
clear pg 
pg.pulses=56;
pg.ctrl='loop pack'; 
burstR={};
plen = 12;
V_rf = vrf;
rtmp = pdload('right');
% Parameters: [pulselength, RF_phase, eps_offset(x,y), eps_RF(x,y), RF_freq, total_evo_time, dt] 
%pg.dict={struct('meas','@meas12','prep','@adprep', 'read','@adread'),'right'};
sepval = rtmp.sep.val;
%offst = sepval-.5*V_rf*(sepval./abs(sepval));
offst = -.5*sepval./abs(sepval)+sepval-.5*V_rf*sepval./abs(sepval);
namepat='rotating_ramseyE_%02d_R'; pg.chan = [3 4]; 
clear RFypi;
rfFreq = rtmp.RFypi.val(1);
% rfFreq = 176; piamp = 3.5; pitime = .03;
% if 1 %1 for no wait before and after the pi pulse
% RFypi.type = '@RFburst';
% RFypi.time = pitime;
% RFypi.val = [rfFreq,piamp*.5*sepval./abs(sepval),offst,pi/2];
% else
% RFypi(1).type = '@wait'; RFypi(1).time = .005; RFypi(1).val = [-6 6];
% RFypi(2).type = '@RFburst';
% RFypi(2).time = pitime;
% RFypi(2).val = [rfFreq,piamp*.5*sepval./abs(sepval),offst,pi/2];
% RFypi(3) = RFypi(1);
% end

pg.dict={struct('prep','@adprep', 'read','@adread','pi','@RFypi'),'right'};

%tau=.2:.2:4; 
tau = linspace(.3,8,20);
for j=tau
    pg.xval = j;
    pg.params = [plen, 0, offst, V_rf*.5*sepval./abs(sepval), rfFreq , j ,0];
    pg.varpar=4*(1:64)'-2*64';
    pg.name=sprintf(namepat,length(burstR)+1);
    try 
         plsupdate(pg);
    catch
         plsdefgrp(pg);
    end
    burstR = [burstR, pg.name];
end

awgadd(burstR);
awgcntrl('on start wait err raw');
if sm_setgradient;
fprintf('dbz = %.2f\n',sm_getgradient);
fprintf('VRF =%.2f\n',vrf);
smrun(fConfSeq2([28:48],struct('nloop',100,'nrep',400, 'opts','pol','datachan','DAQ2')),smnext('ramseyE_rotate_R'));
fprintf('dbz = %.2f \n',sm_getgradient);
else
   fprintf('no lock gradient for vrf = %.2f \n',vrf); 
end
end
sleep


%% rotating frame CPMG2


for eps = linspace(3,4,10)
clear pg 
pg.pulses=62;
pg.ctrl='loop pack'; 
genR={};
plen = 6;
% Parameters: [pulselength, epsPi, pitime, eps, evotime, dt,] 

for j = linspace(.3,3,20);
pg.chan = [3 4]; 
pg.dict={struct('prep','@adprep', 'read','@adread','pi','@RFypi'),'right'};
pg.varpar=(1:64)'-32;%%(1:128)'-64;%128+0*(1:128)';
pg.params = [plen, eps j,0];
pg.xval = j;
pg.name = sprintf('ramseyCPMG_rot_phaseGen__%02d_2013_04_17_R',length(genR)+1);
genR = [genR pg.name];
plsupdate(pg);
end

awgadd(burstR);
awgcntrl('on start wait err raw');
if sm_setgradient;
fprintf('dbz = %.2f\n',sm_getgradient);
fprintf('VRF =%.2f\n',eps);
d=smrun(fConfSeq2([28:48],struct('nloop',100,'nrep',400, 'opts','pol','datachan','DAQ2')),smnext('ramseyCPMG2_rotate_R'));
if any(isnan(d{1}))
    break
end
fprintf('dbz = %.2f \n',sm_getgradient);
else
   fprintf('no lock gradient for vrf = %.2f \n',eps); 
end

end

sleep

%% rotating frame CPMG (3 pulse)

clear pg 
pg.pulses=59;
pg.ctrl='loop pack'; 


plen = 12;
'hi'
for vrf =[3.1 3.3 3.5]
V_rf = vrf;
burstR={};

rtmp = pdload('right');
% Parameters: [pulselength, RF_phase, eps_offset(x,y), eps_RF(x,y), RF_freq, total_evo_time, dt] 

%pg.dict={struct('meas','@meas12','prep','@adprep', 'read','@adread'),'right'};

sepval = rtmp.sep.val;
offst = sepval-.5*V_rf*(sepval./abs(sepval));

namepat='rotating_ramseyCPMG3_%02d_R'; pg.chan = [3 4]; 

clear RFypi;
rfFreq = 176; piamp = 3.5; pitime = .03;
if 1 %1 for no wait before and after the pi pulse
RFypi.type = '@RFburst';
RFypi.time = pitime;
RFypi.val = [rfFreq,piamp*.5*sepval./abs(sepval),offst,pi/2];
else
RFypi(1).type = '@wait'; RFypi(1).time = .005; RFypi(1).val = [-6 6];
RFypi(2).type = '@RFburst';
RFypi(2).time = pitime;
RFypi(2).val = [rfFreq,piamp*.5*sepval./abs(sepval),offst,pi/2];
RFypi(3) = RFypi(1);
end
pg.dict={struct('prep','@adprep', 'read','@adread','pi',RFypi),'right'};

tau = linspace(.3,8,20);
for j=tau
    pg.xval = j;
    pg.params = [plen, 0, offst, V_rf*.5*sepval./abs(sepval), rfFreq , j ,0];
    pg.varpar=4*(1:64)'-2*64';
    pg.name=sprintf(namepat,length(burstR)+1);
    try 
         plsupdate(pg);
    catch
         plsdefgrp(pg);
    end
    burstR = [burstR, pg.name];
end

awgadd(burstR);
awgcntrl('on start wait err raw');
if sm_setgradient;
fprintf('dbz = %.2f\n',sm_getgradient);
fprintf('VRF =%.2f\n',vrf);
d=smrun(fConfSeq2([28:48],struct('nloop',100,'nrep',400, 'opts','pol','datachan','DAQ2')),smnext('ramseyCPMG3_rotate_R'));
if any(isnan(d{1}))
   break 
end
fprintf('dbz = %.2f \n',sm_getgradient);
else
   fprintf('no lock gradient for vrf = %.2f \n',vrf); 
end

end

sleep

%%
%% rotating frame, CPMG, external generator, arbitrary number of pulses

evotimes = linspace(0.3,8,20);
varpar= (1:64)'-32;
eps = 3.45;
plen = 12;
% below here should happen automatically

clear pg 
pg.pulses=66;
pg.ctrl='loop pack'; 

rtmp = pdload('right');
ee = rtmp.markerburst; ee.val(3:7) = [eps, 0, 1, 0, 0];

for npls =[1 2 4 8 10 11 12 20 21 32 33]% [1:4, 6, 8, 10, 12]
    % Parameters: [pulselength, eps, time_of_last_evo]
    genR={};
    for j = evotimes;
        pg.chan = [3 4];
        evol = make_cpmg_evo(ee,rtmp.RFypi,npls,j);
        pg.dict={struct('evo',evol(1:end-1),'prep','@adprep', 'read','@adread'),'right'};
        pg.varpar=(j/(2*npls))+1e-3*varpar;%%(1:128)'-64;%128+0*(1:128)';
        pg.params = [plen, eps 0];
        pg.xval = j;
        pg.name = sprintf('ramseyCPMG%01d_generic_%02d_2013_04_17_R',npls,length(genR)+1);
        genR = [genR pg.name];
        try
            plsupdate(pg);
        catch
            plsdefgrp(pg)
        end
    end
    awgrm(28,'after'); awgclear('unused');
    fprintf('loading CPMG with %01d pulses\n',npls);
    awgadd(genR);
    awgcntrl('on start wait err raw');
    if sm_setgradient;
        fprintf('dbz = %.2f \n',sm_getgradient);
        fprintf('npls =%01d \n',npls);
        d=smrun(fConfSeq2([28:48],struct('nloop',150,'nrep',250, 'opts','pol','datachan','DAQ2')),smnext(sprintf('ramseyCPMG%01d_rotate_R',npls)));
        if any(isnan(d{1}))
            break
        end
        fprintf('dbz = %.2f \n',sm_getgradient);
    else
        fprintf('no lock gradient for npls = %01d \n',npls);
    end
    
end
sleep

%% echo with external generator, phase control

for eps = 3.5%linspace(4,3,10)
clear pg 
pg.pulses=60;
pg.ctrl='loop pack'; 
genR={};

plen = 12;
pitime = .022;
pieps = 3.6452;
rtmp = pdload('right');
bstart = rtmp.wait.time + sum(rtmp.reload.time) +rtmp.adprep.time-.01;
for j = linspace(.3,6,20);
pg.chan = [3 4]; 
pg.dict={struct('prep','@adprep', 'read','@adread'),'right'};
pg.varpar=(1:64)'-32;%%(1:128)'-64;%128+0*(1:128)';
pg.params = [plen, 1e-3*max(pg.varpar), bstart, pitime, pieps, eps, j, 0];
pg.xval = j;
pg.name = sprintf('ramsey_rot_phaseGen_%02d_2013_04_17_R',length(genR)+1);
genR = [genR pg.name];
try
   plsupdate(pg);
catch
    plsdefgrp(pg);
end
end
awgadd(genR);

awgcntrl('on start wait err raw');
if sm_setgradient;
fprintf('dbz = %.2f\n',sm_getgradient);
fprintf('VRF =%.2f\n',eps);
d=smrun(fConfSeq2([28:48],struct('nloop',150,'nrep',100, 'opts','pol','datachan','DAQ2')),smnext('ramseyE_rotate_genPhase_R'));
if any(isnan(d{1}))
   break 
end
fprintf('dbz = %.2f \n',sm_getgradient);
else
   fprintf('no lock gradient for vrf = %.2f \n',eps); 
end

end
sleep

%% eps sweep, chopping ramsey_rotate w dbz
epsvals = [3.9, (3:.1:5)];
for j = epsvals
clear pg 
pg.pulses=[12 63];
pg.ctrl='multi loop pack'; 
genR={};
plen = 4;
pg.chan = [3 4]; 
dbzvarpar = 0*(0:1499)'+22;
rotvarpar = (0:599)';
stdly = .03; eddly = 0;
% Parameters: [pulselength, startdelay, enddelay,  ch1m1, ch1m2, ch2m1, ch2m2, eps, time,] 
%pg.dict={struct('meas','@meas12','prep','@adprep', 'read','@adread'),'right'};
pg.dict={struct('prep','@adprep', 'read','@adread'),'right'};
pg.varpar={dbzvarpar, rotvarpar};%128+0*(1:128)';
pg.params = {[plen 300 NaN 0],[plen, stdly, eddly, 0, 1, 0, 0, j, 0]};
pg.xval = j;
pg.name=sprintf(('rotramsey_chopdbz_%02d_2013_04_29E_R'),length(genR)+1);
genR = [genR, pg.name];
plsupdate(pg);
awgadd(genR);
awgcntrl('on start wait err raw');
if sm_setgradient;
fprintf('dbz = %.2f\n',sm_getgradient);
fprintf('Voffset =%.2f\n',j);
scanname = smnext(sprintf('ramsey_chop_rot_%.2f_R',j));
scan = fConfSeq2_v2([28],struct('nloop',1,'nrep',1000, 'opts','pol','datachan','DAQ2'));
scan.saveloop = [1 20];
d=smrun(scan,scanname);
if any(isnan(d{1}))
   break 
end
fprintf('dbz = %.2f \n',sm_getgradient);
else
   fprintf('no lock gradient for vrf = %.2f \n',j); 
end
end
sleep

