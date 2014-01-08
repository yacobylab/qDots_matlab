%%

% 6us pulses, pol
scanEDSR2_pol = fConfSeq2_v2(30,struct('nloop',50,'nrep',1024,'datachan','DAQ2','opts','pol logdbz'));
scanEDSR2_pol.saveloop = [1 100];
scanEDSR2_pol.loops(2).setchan = 'count'; scanEDSR2_pol.loops(2).npoints = 10; scanEDSR2_pol.loops(2).rng = [1 10];
scanEDSR2_pol.loops(1).setchan = 'RFfreq3';
scanEDSR2_pol.loops(1).rng = [3.5 6.5]*1e9;

%sm_setgradient; sm_getgradient;
smrun(scanEDSR2_pol,smnext('EDSR_eps_6us_pol_R'));
%sleep; 
%pause(120);

% 4us pulses, pol
scanEDSR2_pol = fConfSeq2_v2(31,struct('nloop',75,'nrep',1024,'datachan','DAQ2','opts','pol logdbz'));
scanEDSR2_pol.saveloop = [1 100];
scanEDSR2_pol.loops(2).setchan = 'count'; scanEDSR2_pol.loops(2).npoints = 10; scanEDSR2_pol.loops(2).rng = [1 10];
scanEDSR2_pol.loops(1).setchan = 'RFfreq3';
scanEDSR2_pol.loops(1).rng = [1 4]*1e9;

%sm_setgradient; sm_getgradient;
smrun(scanEDSR2_pol,smnext('EDSR_eps_4us_pol_R'));
%%
sleep; 
pause(600);
%%
%next is 6 us pulses, no pol
scanEDSR2 = fConfSeq2_v2('EDSR_eps_2013_09_15_R',struct('nloop',400,'nrep',512,'datachan','DAQ2','opts','logdbz'));
scanEDSR2.saveloop = [1 100];
scanEDSR2.loops(2).setchan = 'count'; scanEDSR2.loops(2).npoints = 10; scanEDSR2.loops(2).rng = [1 10];
scanEDSR2.loops(1).setchan = 'RFfreq3';
scanEDSR2.loops(1).rng = [4.5 6.5]*1e9;
%%

scanEDSR2 = fConfSeq2_v2('EDSR_eps_2013_09_15_R',struct('nloop',75,'nrep',512,'datachan','DAQ2','opts','logdbz'));
%scanEDSR2 = fConfSeq2_v2(29,struct('nloop',75,'nrep',400,'datachan','DAQ2','opts','logdbz'));
scanEDSR2.saveloop = [1 100];
scanEDSR2.loops(2).setchan = 'count'; scanEDSR2.loops(2).npoints = 6; scanEDSR2.loops(2).rng = [1 6];
scanEDSR2.loops(1).setchan = 'RFfreq3';
scanEDSR2.loops(1).rng = [4 7]*1e9;
d=smrun(scanEDSR2,smnext('EDSR_R'));
sleep
scanEDSR2 = fConfSeq2_v2('EDSR_eps_2013_09_15_R',struct('nloop',75,'nrep',256,'datachan','DAQ2','opts','logdbz'));

scanEDSR2.saveloop = [1 100];
scanEDSR2.loops(2).setchan = 'count'; scanEDSR2.loops(2).npoints = 6; scanEDSR2.loops(2).rng = [1 6];
scanEDSR2.loops(1).setchan = 'RFfreq3';
scanEDSR2.loops(1).rng = [.1 2]*1e9;
d=smrun(scanEDSR2,smnext('EDSR_ST0_R'));
%%
scanEDSR3 = fConfSeq2_v2(31,struct('nloop',400,'nrep',512,'datachan','DAQ2','opts','logdbz'));
scanEDSR3.saveloop = [1 100];
scanEDSR3.loops(2).setchan = 'count'; scanEDSR3.loops(2).npoints = 10; scanEDSR3.loops(2).rng = [1 10];
scanEDSR3.loops(1).setchan = 'RFfreq3';
scanEDSR3.loops(1).rng = [4.5 6.5]*1e9;



%%
pause(350);
scanEDSR2 = fConfSeq2_v2(30,struct('nloop',75,'nrep',400,'datachan','DAQ2','opts','logdbz'));
scanEDSR2.saveloop = [1 100];
scanEDSR2.loops(2).setchan = 'count'; scanEDSR2.loops(2).npoints = 5; scanEDSR2.loops(2).rng = [1 5];
scanEDSR2.loops(1).setchan = 'RFfreq3';
scanEDSR2.loops(1).rng = [3 5]*1e9;
smrun(scanEDSR2,smnext('EDSR_eps_6us_R'));
sleep

%%
fvals = linspace(5.3e9,5.6e9,32);
for j = 1:length(fvals)
   smset('RFfreq3',fvals(j));
   fprintf('setting frequency to %.2dGHz \n',1e-9*fvals(j));
   d=smrun(fConfSeq2(52,struct('nloop',1200,'nrep',500, 'opts','logdbz','datachan','DAQ2')),smnext('EDSR_freq_R'));
    if any(isnan(d{1}))
        sleep
        break
    end
end
sleep
%%
scanEDSR3 = fConfSeq2_v2(31,struct('nloop',200,'nrep',128,'datachan','DAQ2','opts','logdbz'));
scanEDSR3.saveloop = [1 100];
scanEDSR3.loops(2).setchan = 'count'; scanEDSR3.loops(2).npoints = 16; scanEDSR3.loops(2).rng = [1 16];
scanEDSR3.loops(1).setchan = 'RFpow3';
scanEDSR3.loops(1).rng = [-5, 5];

%% sweep evo time (crudely)

clear pg 
pg.pulses=73;
pg.ctrl='loop pack'; 
genR={};
%Parameters: [pulselength, phase, eps_offset(x,y), eps_RF(x,y), RF_freq, burst time] 
plen = 6;
pg.chan = [3 4]; 

%Parameters:  pulselength, time(ns), eps
aa=tic;
for Tevo = [2048]
pg.dict={'right'};
pg.varpar=linspace(1.5,5.5,128)'; 
pg.name=sprintf(('EDSR_eps_2013_09_15_R'),length(genR)+1);
%pg.varpar = [8*(1:128)', ones(1,128)'*5]; pg.name=sprintf(('EDSR_eps_2013_09_16_R'),length(genR)+1);
pg.params = [plen, Tevo,0];
plsupdate(pg);
awgadd(pg.name);
awgcntrl('on start wait err raw');
d=smrun(scanEDSR2,smnext(sprintf('EDSR_eps_%.0d_6us_R',Tevo)));
if any(isnan(d{1}))
    sleep
    break
end
end
fprintf('%.2d\n',toc(aa));
sleep

%% fast passage

scanEDSR2 = fConfSeq2_v2('arp_adeps_2013_10_11_R',struct('nloop',40,'nrep',256,'datachan','DAQ2','opts','logdbz'));
%scanEDSR2 = fConfSeq2_v2('EDSR_eps_2013_09_15_R',struct('nloop',400,'nrep',256,'datachan','DAQ2','opts','logdbz'));

scanEDSR2.saveloop = [1 100];
scanEDSR2.loops(2).setchan = 'count'; scanEDSR2.loops(2).npoints = 10; scanEDSR2.loops(2).rng = [1 10];
scanEDSR2.loops(1).setchan = 'RFfreq3';
scanEDSR2.loops(1).rng = [4.75 6.25]*1e9;

%% Coherent beam splitting
%stpscan=fConfSeq2([30:61],struct('nloop',500, 'nrep',100,'opts','','datachan','DAQ2'));
stppolscan=fConfSeq2(28,struct('nloop',1000, 'nrep',500,'opts','','datachan','DAQ2'));
%stpscan=fConfSeq2([30:61],struct('nloop',500, 'nrep',100,'opts','','datachan','DAQ2'));


%% STp spectroscopy
stpscan=fConfSeq2(61,struct('nloop',256, 'nrep',128,'opts','','datachan','DAQ2'));
stpscan.saveloop = [2 1];
stpscan.loops(1).setchan = 'RFfreq3';
stpscan.loops(1).rng = [10e6 6e9];
stpscan.loops(2).setchan = 'count';
stpscan.loops(2).npoints =50;
stpscan.loops(2).rng = [1 stpscan.loops(2).npoints];

testscan=fConfSeq2(61,struct('nloop',2000, 'nrep',64,'opts','','datachan','DAQ2'));
