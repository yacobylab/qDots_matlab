%% Now do a tomography calibration

clear opts;
opts.lastgroup = 27;
opts.ffunc_i=str2func('@(p,x) p(2)*log(p(1)./(x-p(3)))');

compensate = 1;
awgrm(opts.lastgroup,'after');
awgclear('unused');
stag=0; %do staggered read out. 
for j =1;%:1

if j %right side
    opts.chan = [3 4];
    opts.side = 'right';
    opts.skin=0;
    opts.rc=0;
    opts.ramp = 0;
    opts.datachan = 'DAQ2';
    opts.fp=[5.0986e+03 1.0620 0];%[1.3e+03 0.65 0];
    opts.oside = 'left';
    namepat= 'tomocal_2012_12_08_R_rot_%d_%s';
else %left side
    opts.chan = [2 1];
    opts.side = 'left';
    opts.oside = 'right';
    opts.skin=0;
    opts.rc=0;
    opts.ramp = 0;
    opts.datachan = 'DAQ1';
    opts.fp=[6.0384e+03 1.1860 0];%[7.8632e+03 0.6767 0];%[3.9592e+03 0.9377 0];
    namepat= 'tomocal_2012_12_08_L_rot_%d_%s';
end

%awgrm(opts.lastgroup,'after');
%awgclear('unused');
clear tcal;
slowfudge=.1;
fastfudge=0;
sepper = 30;
tcal(1).prep='@null'; tcal(1).per =sepper; tcal(1).repeat=0;
tcal(2).prep='@dbzprep'; tcal(2).per=sepper/4; tcal(2).repeat=1;
tcal(3).prep='@adprep'; tcal(3).per=sepper/6; tcal(3).repeat=1;
tcal(4).prep='@dbzprep'; tcal(4).per=15; tcal(4).repeat=0;
tcal(5).prep='@dbzpi'; tcal(5).per=25; tcal(5).repeat=0;
tcal(6).prep='@null'; tcal(6).per=32*5; tcal(6).repeat=0;
tomos={'ST','UD','Y'};
tcalg={};
offset=0;
clear pg;

if stag
      rd = ['stag', opts.side(1)];
      opts.side2 = {rd, opts.side};
else
    opts.side2=opts.side;
end
  
for i=1:length(tcal)
    for j=1:length(tomos)
        pg.ctrl='loop pack';
        pg.chan = opts.chan;
        pg.dict={struct('prep',tcal(i).prep,'tomo',['@', tomos{j}, 'read']), opts.side2};
        pg.name=sprintf(namepat,i,tomos{j});
        pg.pulses=34;
        pg.params=[0 0];
        pg.trafofn(1).func=@skineffect_trafofn;
        pg.trafofn(1).args=opts.skin;
        pg.trafofn(2).func=@rc_trafofn;
        pg.trafofn(2).args=opts.rc;
        pg.trafofn(3).func=@rampify_trafofn;
        pg.trafofn(3).args=opts.ramp;
        
        if tcal(i).repeat
            pg.varpar=[(1+offset):(16+offset) (1+offset):(16+offset)]' -1;
        else
            pg.varpar=((1+offset):(32+offset))' -1;
        end
        if tcal(i).per < 32
            fudge=fastfudge;
        else
            fudge=slowfudge;
        end

        pg.params(1)=opts.ffunc_i(opts.fp,1e3/tcal(i).per) + fudge;
        pg.params;
        if ~isreal(pg.params(1))
            pg.params(1)=out.jcal.ffunc_i(out.jcal.fp .* [1 1 0],1e3/tcal(i).per) + fudge;
        end
        try
            plsupdate(pg);
        catch e;
            fprintf('Error adding group: %s: %s\n',e.identifier,e.message);
            plsdefgrp(pg);
        end
        tcalg=[tcalg pg.name];
    end
end

if compensate
    clear pg;
    %namepat = 'tomocalrot_%d_R_zero_L_';
    comptcalg = {};
    load compmatrix;
    pg.matrix = compmatrix;
    compJgrp = {};
    pg.chan=[1 2 3 4];
   zgp = sprintf('zeros_10_32_%s',upper(opts.oside(1)));
    for j=1:length(tcalg)
        pg.pulses.groups={tcalg{j},zgp};
        pg.name=[tcalg{j},sprintf('_zero_%s',upper(opts.oside(1)))];
        comptcalg = [comptcalg, pg.name];
        pg.ctrl='grp loop pack';
        pg.pulseind(2,:)=[1:32];
        pg.pulseind(1,:)=[1:32];
        try
          plsupdate(pg);
        catch
           plsdefgrp(pg); 
        end
        
    end
    fprintf('loading tomocal pulses\n')
   awgadd(comptcalg);
else
    fprintf('loading tomocal pulses\n')
    awgadd(tcalg);
end

adddbz(upper(opts.side(1)));

end
awgcntrl('on start wait err raw');
% analyze with ana_Tomocal_2

%smrun(fConfSeq2([45 27:44],struct('nloop',200,'nrep',40,'opts','pol','datachan','DAQ1')),smnext('tomocal_Rot_L_zero_R'))
%smrun(fConfSeq2([64 46:63],struct('nloop',200,'nrep',40,'opts','pol','datachan','DAQ2')),smnext('tomocal_Rot_R_zero_L'))
