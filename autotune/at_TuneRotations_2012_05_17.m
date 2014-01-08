function [ out ] = at_TuneRotations(opts, out)
% function [ out ] = at_TuneRotations(opts, out)
%  Tune up a set of rotations.  out has a record of the "thought" process
%  of the rotation tuner, and can be passed back in to have at_TuneRotations
%  add some extra rotations to a tuned set.
% at_TuneRotations uses tomography based tunings of rotations.
%
% WARNING: This function adds and deletes groups, reprograms the AWG,
%  restarts the AWG, calls autotune, locks the gradient, and recenters the
%  dot using atcenter
%
% options:
% opts.opts: string containing
%   ana -- analyze only, don't take data.
% opts.lastgroup: TuneRotations will delete groups after this.  Default 23
% opts.side: default "right"
% opts.chan: default [3 4]
% opts.tests: string containing
%      t1 -- check t1
% tomocal -- calibrates tomography.  Implies t1
%

    global tunedata;
    global fbdata;
    global awgdata;
    
    if ~exist('opts','var')
        opts=struct();
    else
        if isfield(opts,'log')
            error('Did you pass an old output to the options?');
        end        
        if ischar(opts)
            opts=struct('opts',opts);
        end
    end
    try

    opts=def(opts,'opts','');  % ana
    ana = isopt(opts,'ana');

    opts=def(opts,'tests','tunes');
    opts=def(opts,'side','right');
    opts=def(opts,'chan',[3 4]);       % fixme; get from tunedata
    opts=def(opts,'datachan','DAQ2');  % fixme; get from tunedata
    opts=def(opts,'lastgroup',23);
    opts=def(opts,'skin',0); % Skin effect @ 1 Ghz in dB
    opts=requires(opts,'all',{'tests','tunes','proccal_load'});
    opts=requires(opts,'tests',{'proctest_j_plus','proctest_j_minus','proctest_j_pi'});
    %opts=requires(opts,'tunes',{'proccal_j_pi','proccal_hjh','proccal_hjjh','proccal_j_plus','proccal_j_minus','proccal_jjh','proccal_hjmh'});
    opts=requires(opts,'tunes',{'proccal_j_pi','proccal_hjh','proccal_hjjh','proccal_j_plus','proccal_j_minus','proccal_jjh'});
    if ~ana
        opts=requires(opts,'proccal_jjh',{'proccal_j_pi','proccal_hjh','proccal_hjmh'});
        opts=requires(opts,'proccal_hjh','proccal_j_plus');
        opts=requires(opts,'proccal_hjjh','proccal_j_pi');
        %opts=requires(opts,'proccal_hjmh','proccal_j_minus');
        opts=requires(opts,'proctest_j_plus','proccal_j_plus');
        opts=requires(opts,'proccal_j_plus','tomocal');
        opts=requires(opts,'proctest_j_pi','proccal_j_plus');
        opts=requires(opts,'proccal_j_pi','tomocal');
        opts=requires(opts,'proctest_j_minus','proccal_j_minus');
        opts=requires(opts,'proccal_j_minus','tomocal');
        opts=requires(opts,'proccal_b_plus','tomocal');
        opts=requires(opts,'proccal_load','tomocal');
        opts=requires(opts,'tomocal',{'t1','jcal'});
        opts=requires(opts,'jcal','t1');
    end    
    
    if ~exist('out','var')
        out.log=sprintf('Started a new at_TuneRotations run at %s\n',datestr(now));
        out.fname=smnext('attr_summary');
        out.figs=[];
    else
        out.log=[out.log sprintf('Continued at_TuneRotations run at %s\n',datestr(now))];
    end
    out.log=[out.log sprintf('Tests to run: %s\n%g dB skin comp\n',opts.tests,opts.skin)];
    
    % If ana is set, don't run any scans or change any gates.
    ana = isopt(opts,'ana');
    if ana
        anastr=' ana';
    else
        anastr='';
    end
    
    dict=pdload(opts.side);

    
    %% Start out by checking T1
    if istest(opts,'t1') && (ana || ~isdone(out,'t1'))
        if ana
            autotune(['t1' anastr],out.t1.run);
        else
            PrepareMeasurement(opts);
            autotune('t1');
            if ~VerifyLock(opts)
                attr_log(out, '  -> Gradient came unlocked during test!  aborting\n');
                %out.log=[out.log sprintf('  -> Gradient came unlocked during test!  aborting\n')];
                fprintf('Gradient came unlocked during T1 test\n');
                out=cleanup(out,opts);
                return;
            end
            sleep;
            out.t1.run = length(tunedata.runs);
            out.t1.subrun = length(tunedata.runs(end).t1);
        end
        out.t1.vis = tunedata.runs(out.t1.run).t1(out.t1.subrun).vis;
        out.t1.t1  = tunedata.runs(out.t1.run).t1(out.t1.subrun).t1;
        
        if ~ana
            attr_log(out,'  -> Tunerun %d.%d: measured T1=%g, vis=%g\n',[out.t1.run,out.t1.subrun, out.t1.t1, out.t1.vis']);
            %out.log = [out.log sprintf('  -> Tunerun %d.%d: measured T1=%g, vis=%g\n',out.t1.run,out.t1.subrun, out.t1.t1, out.t1.vis')];
            out.t1.done=1;
        end
        out.figs=[out.figs [79 80 81]];
    end
    
    %% Now do an exchange calibration
    if istest(opts,'jcal') && (ana || ~isdone(out,'jcal'))
        %% Load the groups
        if ~ana 
            awgrm(opts.lastgroup,'after');
            clear pg;
            pg.pulses=21804;
            namepat='at_tunerot_jcal_%02d'; pg.chan = opts.chan; pg.dict = opts.side;
            pg.ctrl='loop pack';
            jgrp={};
            pg.dict={struct('prep',struct('type','@dbzprep'),'read',struct('type','@dbzread')),pg.dict};
            for eps=.7:.1:2.5
                pg.xval = eps;
                pg.params=[eps 0];
                pg.varpar=[0:1:63]';
                pg.name=sprintf(namepat,length(jgrp)+1);
                try
                    plsupdate(pg);
                catch
                    plsdefgrp(pg);
                end
                jgrp=[jgrp pg.name];
            end
            awgadd(jgrp);
            adddbz(upper(opts.side(1)));
            
            PrepareMeasurement(opts);
            
            out.jcal.fname=smnext(sprintf('attr_Ramsey_%s',upper(opts.side(1))));
            d=smrun(fConfSeq2([length(awgdata.pulsegroups),(opts.lastgroup+1):(length(awgdata.pulsegroups)-1)], ...
                struct('nloop',400,'nrep',10,'opts','pol','datachan',opts.datachan)),out.jcal.fname);  
            out.jcal.fname=['sm_' out.jcal.fname];
            attr_log(out,'  -> Took a Ramsey dataset (%s)\n',out.jcal.fname);
            %out.log = [out.log sprintf('  -> Took a Ramsey dataset (%s)\n',out.jcal.fname)];
            if any(isnan(d{1}(:)))
                attr_log(out,'    Scan aborted by user\n');
                %out.log=[out.log sprintf('    Scan aborted by user\n')];
                out=cleanup(out,opts);
                return
            end
            if ~VerifyLock(opts)
                attr_log(out,'  -> Gradient came unlocked during test!  aborting\n');
                %out.log=[out.log sprintf('  -> Gradient came unlocked during test!  aborting\n')];
                %fprintf('Gradient came unlocked during Ramsey test\n');
                out=cleanup(out,opts);
                return;
            end
            
        end
        [figs fits data]=ana_echo(out.jcal.fname,struct('opts','ramsey freq afitdecay','rng',[6 inf],'fb',300));
        gooddata=find(fits.freq(2,:) > 200,1,'last'):size(fits.freq,2);   % These are in MHz.
        out.jcal.fits=fits;
        out.jcal.data=data;
        out.jcal.freqs=fits.freq(:,gooddata);
        out.jcal.freqs(2,:)=out.jcal.freqs(2,:)*1e-3;
        out.jcal.ffunc=str2func('@(p,x) p(1)*exp(-x/p(2))+p(3)');
        out.jcal.ffunc_i=str2func('@(p,x) p(2)*log(p(1)./(x-p(3)))');
        out.jcal.done=1;
        subplot(2,2,4);
        out.jcal.fp=fitwrap('samefig plinit plfit',out.jcal.freqs(1,:),out.jcal.freqs(2,:),[.300,.3,0],out.jcal.ffunc);
        xlabel('\epsilon');
        ylabel('F (Ghz)');
        out.figs = [figs out.figs];
    end
    
    %% Now do a tomography calibration
    if istest(opts,'tomocal') && (ana || ~isdone(out,'tomocal'))
        if ~ana
            awgrm(opts.lastgroup,'after');
            clear tcal;
            slowfudge=[.1];
            fastfudge=[0];
            tcal(1).prep='@null'; tcal(1).per =32; tcal(1).repeat=0;
            tcal(2).prep='@dbzprep'; tcal(2).per=32/4; tcal(2).repeat=1;
            tcal(3).prep='@adprep'; tcal(3).per=32/6; tcal(3).repeat=1;
            tcal(4).prep='@dbzprep'; tcal(4).per=15; tcal(4).repeat=0;
            tcal(5).prep='@dbzpi'; tcal(5).per=25; tcal(5).repeat=0;
            tcal(6).prep='@null'; tcal(6).per=32*5; tcal(6).repeat=0;
            tomos={'ST','UD','Y'};
            tcalg={};
 
            clear pg;            
            for i=1:length(tcal)
                for j=1:length(tomos)
                    pg.ctrl='loop pack';
                    pg.chan = opts.chan;
                    pg.dict={struct('prep',tcal(i).prep,'tomo',['@', tomos{j}, 'read']), opts.side};
                    pg.name=sprintf('at_tomocal_rot_%d_%s',i,tomos{j});
                    pg.pulses=21796;                 
                    pg.params=[0 0];
                    pg.trafofn.func=@skineffect_trafofn;
                    pg.trafofn.args=opts.skin;
                    
                    if tcal(i).repeat
                        pg.varpar=[1:16 1:16]' -1;
                    else
                        pg.varpar=[1:32]' -1;
                    end
                    if tcal(i).per < 32
                        fudge=fastfudge;
                    else
                        fudge=slowfudge;
                    end
                    pg.params(1)=out.jcal.ffunc_i(out.jcal.fp,1/tcal(i).per) + fudge;
                    if ~isreal(pg.params(1))
                        pg.params(1)=out.jcal.ffunc_i(out.jcal.fp .* [1 1 0],1/tcal(i).per) + fudge;
                    end
                    try
                        plsupdate(pg);
                    catch
                        plsdefgrp(pg);
                    end
                    tcalg=[tcalg pg.name];
                end
            end
            awgadd(tcalg);        
            adddbz(upper(opts.side(1)));            
            PrepareMeasurement(opts);
            out.tomocal.fname=smnext(sprintf('attr_TomoCal_rot_%s',upper(opts.side(1))));
            d=smrun(fConfSeq2([length(awgdata.pulsegroups),(opts.lastgroup+1):(length(awgdata.pulsegroups)-1)], ...
                struct('nloop',400,'nrep',10,'opts','pol','datachan',opts.datachan)),out.tomocal.fname);  
            out.tomocal.fname=['sm_' out.tomocal.fname];
            attr_log(out,'  -> Took a Tomography Calibration dataset (%s)\n',out.tomocal.fname);
            %out.log = [out.log sprintf('  -> Took a Tomography Calibration dataset (%s)\n',out.tomocal.fname)];
            if any(isnan(d{1}(:)))
                attr_log(out,'    Scan aborted by user\n');
                %out.log=[out.log sprintf('    Scan aborted by user\n')];
                out=cleanup(out,opts);
                return
            end
            if ~VerifyLock(opts)
                attr_log(out,'  -> Gradient came unlocked during test!  aborting\n');
                %out.log=[out.log sprintf('  -> Gradient came unlocked during test!  aborting\n')];
                fprintf('Gradient came unlocked during Tomocal test\n');
                out=cleanup(out,opts);
                return;
            end            
        end
        out.tomocal.cal=ana_tomocal_2(out.tomocal.fname);
        out.figs=[out.figs out.tomocal.cal.figs];
        out.tomocal.done=1;
    end
    
    

    
    %% Now tune up %J-rotations with process tomo
    if istest(opts,'proccal_load') && (ana || ~isdone(out,'proccal_load'))
        % params are pulse length, epsilon ,pulse time
        varpar=linspace(0,.1,65)';
        
        popts=struct('name','proccal_load','process',struct('type','@reload'),'params',[nan nan .1],'varpar',varpar,'pulsenum',21857);
        out=RunProcessTomo(out, opts, popts);
        
        dummy  = ana_processTomo_c(out.proccal_load.fname,struct('tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
        target = [1 0 0 0; 0 0 0 0; 0 0 0 0; 1 0 0 0]; % Load singlet.
        out.proccal_load.pdata = ana_process_matrix(dummy, {'opts','eigenbasis', 'tgt',target});
        out.proccal_load.done=1;
    end
    
    name='proccal_j_plus';
    target = [1 0 0 0; 0 0 -1 0; 0 1 0 0; 0 0 0 1]; %pi/2 around J
    if istest(opts,name) && (ana || ~isdone(out,name))
        % params are pulse length, epsilon ,pulse time
        proctime=2e-3;
        eps = [0, linspace(out.jcal.ffunc_i(out.jcal.fp,70e-3),...
            out.jcal.ffunc_i(out.jcal.fp,300e-3), 64)]';
        varpar=[eps, [0 ; repmat(proctime, length(eps)-1,1)]];
        
        popts=struct('name',name,'process',struct('type','@exch'),'params',[8 1.15 proctime],'varpar',varpar);
        
        out=RunProcessTomo(out, opts, popts);
        
        % Extract the process matrix
        pmat = ana_processTomo_c(out.(name).fname,struct('tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
       
        % Analyze the process into a rotation and dephasing
        pdesc=ana_process_matrix(pmat, {'opts','eigenbasis', 'tgt',target});
        [maxfid indmaxfid]=max([pdesc.pdata.fid]);        
        out.(name).data = pdesc.pdata(indmaxfid);
        out.(name).data.eps = eps(indmaxfid);
        out.(name).data.fid = maxfid;
        out.(name).data.proctime = proctime;
        out.(name).done=1;
        out.(name).target=target;
    end
    tname='proctest_j_plus';
    if istest(opts,tname) && (ana || ~isdone(out,tname))
        out=RunProcessTest(out, name, tname);  
    end
    
    
    name='proccal_j_pi';   
    target = [1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1]; %pi around J
    if istest(opts,name) && (ana || ~isdone(out,name))
        % params are pulse length, epsilon ,pulse time
        proctime=3e-3;
        eps = [0, linspace(out.jcal.ffunc_i(out.jcal.fp,80e-3),...
            out.jcal.ffunc_i(out.jcal.fp,320e-3), 64)]';
        varpar=[eps, [0 ; repmat(proctime, length(eps)-1,1)]];
        
        popts=struct('name',name,'process',struct('type','@exch'),'params',[8 1.15 proctime],'varpar',varpar);
        
        out=RunProcessTomo(out, opts, popts);
        
        % Extract the process matrix
        pmat = ana_processTomo_c(out.proccal_j_plus.fname,struct('tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
        % Analyze the process into a rotation and dephasing
        pdesc=ana_process_matrix(pmat, {'opts','eigenbasis', 'tgt',target});
        [maxfid indmaxfid]=max([pdesc.pdata.fid]);        
        out.(name).data = pdesc.pdata(indmaxfid);
        out.(name).data.eps = eps(indmaxfid);
        out.(name).data.fid = maxfid;
        out.(name).data.proctime = proctime;
        out.(name).done=1;
        out.(name).target=target;
    end
    tname='proctest_j_pi';
    if istest(opts,tname) && (ana || ~isdone(out,tname))
        out=RunProcessTest(out, name, tname);  
    end
    
    name='proccal_j_minus';
    target = [1 0 0 0; 0 0 1 0; 0 -1 0 0; 0 0 0 1]; %-pi/2 around J
    if istest(opts,name) && (ana || ~isdone(out,name))
        % params are pulse length, epsilon ,pulse time
        proctime=4e-3;
        eps = [0, linspace(out.jcal.ffunc_i(out.jcal.fp,140e-3),...
            out.jcal.ffunc_i(out.jcal.fp,600e-3), 64)]';
        varpar=[eps, [0 ; repmat(proctime, length(eps)-1,1)]];
        
        popts=struct('name',name,'process',struct('type','@exch'),'params',[8 1.15 proctime],'varpar',varpar);
        out=RunProcessTomo(out, opts, popts);
        
        pmat  = ana_processTomo_c(out.(name).fname,struct('tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
        pdesc=ana_process_matrix(pmat, {'opts','eigenbasis', 'tgt',target});
        [maxfid indmaxfid]=max([pdesc.pdata.fid]);
        out.(name).data = pdesc.pdata(indmaxfid);
        out.(name).data.eps = eps(indmaxfid);
        out.(name).done=1;
        out.(name).data.proctime=proctime;
        out.(name).data.fid = maxfid;
        out.(name).target=target;        
    end
    tname='proctest_j_minus';
    if istest(opts,tname) && (ana || ~isdone(out,tname))
        out=RunProcessTest(out, name, tname);  
    end
           
    name='proccal_hjh';
    if istest(opts,name) && (ana || ~isdone(out,name))        
        % params are pulse length, epsilon ,pulse time
        numjs = 16; %number of Js to try (loop iterations)]
        numjps= 16;
        jrng = .3;
        jprng = .3;
        jtime = out.proccal_j_plus.data.proctime;
        jptime = 1e-3*round(1e3*dict.dbzpi.time/sqrt(2));

        jval  = linspace(out.proccal_j_plus.data.eps-.75*jrng, out.proccal_j_plus.data.eps+.25*jrng,numjs);
        jpval = linspace(out.jcal.ffunc_i(out.jcal.fp,35e-3)-.5*jprng, out.jcal.ffunc_i(out.jcal.fp,35e-3)+.5*jprng,numjps);
        [j jp]=meshgrid(jval,jpval);
        epsj   = [0; j(:)];        
        epsjp  = [0; jp(:)];        
        jtime  = [0;  repmat(jtime,length(epsj)-1,1)];
        jptime = [0; repmat(jptime,length(epsj)-1,1)];
        varpar = [ epsjp, jptime, epsj, jtime, epsjp, jptime ];
            
        popts=struct('name',name,'nrep', 15,...
                     'process',struct('type','@exch'),...
                     'process2',struct('type','@exch'), ...
                     'process3',struct('type','@exch'), ...
                     'params',[5 1.15 1e-3 1 1e-3 1.15 1e-3],'varpar',varpar,'pulsenum',21859);
        out=RunProcessTomo(out, opts, popts); 
        out=attr_ana_jprime(out,name,[0 0 1]);
        
        out.(name).jeps=interp1(jval,out.(name).bestx);
        out.(name).heps=interp1(jpval,out.(name).besty);
        out.(name).done=1;
    end      
   
    name='proccal_hjmh';
    if istest(opts,name) && (ana || ~isdone(out,name))        
        % params are pulse length, epsilon ,pulse time
        numjs = 16; %number of Js to try (loop iterations)]
        numjps= 16;
        jrng = .4;
        jprng = .4;
        jtime = out.proccal_j_minus.data.proctime;
        jptime = 1e-3*round(1e3*dict.dbzpi.time/sqrt(2));

        jval  = linspace(out.proccal_j_minus.data.eps-.5*jrng, out.proccal_j_minus.data.eps+.5*jrng,numjs);
        jpval = linspace(out.jcal.ffunc_i(out.jcal.fp,35e-3)-.5*jprng, out.jcal.ffunc_i(out.jcal.fp,35e-3)+.5*jprng,numjps);
        [j jp]=meshgrid(jval,jpval);
        epsj   = [0; j(:)];        
        epsjp  = [0; jp(:)];        
        jtime  = [0;  repmat(jtime,length(epsj)-1,1)];
        jptime = [0; repmat(jptime,length(epsj)-1,1)];
        varpar = [ epsjp, jptime, epsj, jtime, epsjp, jptime ];
            
        popts=struct('name',name,'nrep', 15,...
                     'process',struct('type','@exch'),...
                     'process2',struct('type','@exch'), ...
                     'process3',struct('type','@exch'), ...
                     'params',[5 1.15 1e-3 1 1e-3 1.15 1e-3],'varpar',varpar,'pulsenum',21859);
        out=RunProcessTomo(out, opts, popts); 
        out=attr_ana_jprime(out,name,[0 0 -1]);
        out.(name).jeps=interp1(jval,out.(name).bestx);
        out.(name).heps=interp1(jpval,out.(name).besty);
        out.(name).done=1;
    end
    
    name='proccal_hjjh';
    if istest(opts,name) && (ana || ~isdone(out,name))        
        % params are pulse length, epsilon ,pulse time
        numjs = 16; %number of Js to try (loop iterations)]
        numjps= 16;
        jrng = .3;
        jprng = .3;
        jtime = out.proccal_j_pi.data.proctime;
        % For stability it's beneficial to make J' too close to Z axis.
        jptime = 2e-3+1e-3*round(1e3*dict.dbzpi.time/sqrt(2));

        jval  = linspace(out.proccal_j_pi.data.eps-.5*jrng, out.proccal_j_pi.data.eps+.5*jrng,numjs);
        jpval = linspace(out.jcal.ffunc_i(out.jcal.fp,35e-3)-.5*jprng, out.jcal.ffunc_i(out.jcal.fp,35e-3)+.5*jprng,numjps);
        [j jp]=meshgrid(jval,jpval);
        epsj   = [0; j(:)];        
        epsjp  = [0; jp(:)];        
        jtime  = [0;  repmat(jtime,length(epsj)-1,1)];
        jptime = [0; repmat(jptime,length(epsj)-1,1)];
        varpar = [ epsjp, jptime, epsj, jtime, epsjp, jptime ];
            
        popts=struct('name',name,'nrep', 15,...
                     'process',struct('type','@exch'),...
                     'process2',struct('type','@exch'), ...
                     'process3',struct('type','@exch'), ...
                     'params',[5 1.15 1e-3 1 1e-3 1.15 1e-3],'varpar',varpar,'pulsenum',21859);
        out=RunProcessTomo(out, opts, popts); 
 
        target=eye(4);
        target(2:end,2:end)=rmat(out.proccal_hjh.tgtaxis+out.proccal_hjmh.tgtaxis,pi);
   
        pmat  = ana_processTomo_c(out.(name).fname,struct('tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
        pdesc=ana_process_matrix(pmat, {'opts','eigenbasis', 'tgt',target});
        [maxfid indmaxfid]=max([pdesc.pdata.fid]);
        out.(name).data = pdesc.pdata(indmaxfid);
        out.(name).data.epsj = epsj(indmaxfid);
        out.(name).data.epsjp = epsjp(indmaxfid);
        out.(name).jtime = jtime;
        out.(name).jptime = jptime; 
        out.(name).data.fid = maxfid;
        out.(name).target=target;        
        out.(name).done=1;
    end    
    
    name='proccal_jjh';  % j-pi followed by H- should give Y pi/2
    if istest(opts,name) && (ana || ~isdone(out,name))        
        % params are pulse length, epsilon ,pulse time
        numjs = 16; %number of Js to try (loop iterations)]
        numjps= 16;
        jrng = .3;
        jprng = .3;
        jtime = out.proccal_j_pi.data.proctime;
        % For stability it's beneficial to make J' too close to Z axis.
        jptime = 2e-3+1e-3*round(1e3*dict.dbzpi.time/sqrt(2));

        jval  = linspace(out.proccal_j_pi.data.eps-.5*jrng, out.proccal_j_pi.data.eps+.5*jrng,numjs);
        jpval = linspace(out.jcal.ffunc_i(out.jcal.fp,35e-3)-.5*jprng, out.jcal.ffunc_i(out.jcal.fp,35e-3)+.5*jprng,numjps);
        [j jp]=meshgrid(jval,jpval);
        epsj   = [0; j(:)];        
        epsjp  = [0; jp(:)];        
        jtime  = [0;  repmat(jtime,length(epsj)-1,1)];
        jptime = [0; repmat(jptime,length(epsj)-1,1)];
        varpar = [ epsj, jtime, epsjp, jptime ];
            
        popts=struct('name',name,'nrep', 15,...
                     'process',struct('type','@exch'),...
                     'process2',struct('type','@exch'), ...                     
                     'params',[5 1.15 1e-3 1 1e-3],'varpar',varpar,'pulsenum',21858);
        out=RunProcessTomo(out, opts, popts); 
        target=eye(4);
        target(2:end,2:end)=rmat(cross([0 0 1],out.proccal_hjh.tgtaxis+out.proccal_hjmh.tgtaxis),pi/2);
        
        pmat  = ana_processTomo_c(out.(name).fname,struct('tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
        pdesc=ana_process_matrix(pmat, {'opts','eigenbasis', 'tgt',target});
        [maxfid indmaxfid]=max([pdesc.pdata.fid]);
        out.(name).data = pdesc.pdata(indmaxfid);
        out.(name).data.epsj = epsj(indmaxfid);
        out.(name).data.epsjp = epsjp(indmaxfid);
        out.(name).jtime = jtime;
        out.(name).jptime = jptime; 
        out.(name).data.fid = maxfid;
        out.(name).target=target;        
        out.(name).done=1;
    end    
       
    name='proccal_hjj';  % j-pi followed by H- should give Y pi/2
    if istest(opts,name) && (ana || ~isdone(out,name))        
        % params are pulse length, epsilon ,pulse time
        numjs = 16; %number of Js to try (loop iterations)]
        numjps= 16;
        jrng = .3;
        jprng = .3;
        jtime = out.proccal_j_pi.data.proctime;
        % For stability it's beneficial to make J' too close to Z axis.
        jptime = 2e-3+1e-3*round(1e3*dict.dbzpi.time/sqrt(2));

        jval  = linspace(out.proccal_j_pi.data.eps-.5*jrng, out.proccal_j_pi.data.eps+.5*jrng,numjs);
        jpval = linspace(out.jcal.ffunc_i(out.jcal.fp,35e-3)-.5*jprng, out.jcal.ffunc_i(out.jcal.fp,35e-3)+.5*jprng,numjps);
        [j jp]=meshgrid(jval,jpval);
        epsj   = [0; j(:)];        
        epsjp  = [0; jp(:)];        
        jtime  = [0;  repmat(jtime,length(epsj)-1,1)];
        jptime = [0; repmat(jptime,length(epsj)-1,1)];
        varpar = [ epsjp, jptime, epsj, jtime ];
            
        popts=struct('name',name,'nrep', 15,...
                     'process',struct('type','@exch'),...
                     'process2',struct('type','@exch'), ...                     
                     'params',[5 1.15 1e-3 1 1e-3],'varpar',varpar,'pulsenum',21858);
        out=RunProcessTomo(out, opts, popts); 
        target=eye(4);
        target(2:end,2:end)=rmat(cross([0 0 1],out.proccal_hjh.tgtaxis+out.proccal_hjmh.tgtaxis),-pi/2);
        
        pmat  = ana_processTomo_c(out.(name).fname,struct('tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
        pdesc=ana_process_matrix(pmat, {'opts','eigenbasis', 'tgt',target});
        [maxfid indmaxfid]=max([pdesc.pdata.fid]);
        out.(name).data = pdesc.pdata(indmaxfid);
        out.(name).data.epsj = epsj(indmaxfid);
        out.(name).data.epsjp = epsjp(indmaxfid);
        out.(name).jtime = jtime;
        out.(name).jptime = jptime; 
        out.(name).data.fid = maxfid;
        out.(name).target=target;        
        out.(name).done=1;
    end    
    
    
catch e
    fprintf('Error running tunerotation.  See out.lasterror\n');  
    fprintf('%s: %s\n',e.identifier,e.message);    
    out.lasterror=e;
    out=cleanup(out,opts);
    return;
    end
out=cleanup(out,opts);
end

function d=findzeroc(data) % find zero crossings
  d=contourc(data,[0 0]);
  ind=1;
  if isempty(d)
      return;
  end
  while 1
      dind=d(2,ind);
      d(:,ind)=nan;
      ind=ind+dind+1;
      if ind > size(d,2)
          break;
      end
  end
end

%% Analyze two-angle rotations
function out = attr_ana_jprime(out,name,paxis)
% function out = attr_ana_jprime(out,name,paxis)
%  paxis is axis to be perpendicular to.
  file=out.(name).fname;
  tomofile=[out.tomocal.fname '_mat.mat'];
  
  process=ana_processTomo_c(file,struct('tomocal',tomofile,'opts',''));
  matrix=ana_process_matrix(process,{'opts','eigenbasis'});
  
  results=matrix.pdata(2:end);
  for i=1:length(results)
      angle(i)=results(i).eigs.angle;
      axis(i,:)=results(i).eigs.axis;
      dp(i)=dot(results(i).eigs.axis,paxis);  % Fixme; compare to j-rotation?
  end
  s=sqrt(length(angle));
  angle=reshape(angle,s,s);
  dp=reshape(dp,s,s);
  %%
  angle=angle-pi/2;
  cang=findzeroc(angle);
  cdp =findzeroc(dp);

  figure(221); clf;
  subplot(121);
  imagesc(angle*180/pi); hold on; title('Rotation Error (degrees)')
  plot(cang(1,:),cang(2,:),'m-');
  cmap_center; colormap(cmap_str('kr.gk')); colorbar;
  subplot(122);
  imagesc(dp); hold on; colorbar; title('Axis Error');
  cmap_center; colormap(cmap_str('kr.gk')); colorbar;
  plot(cdp(1,:),cdp(2,:),'m-');
  

  [xp,yp]=intersections(cang(1,:),cang(2,:),cdp(1,:),cdp(2,:));
  if isempty(xp)
      error('No optimal rotation');
  end
  for i=1:length(xp)
      x=xp(i);
      y=yp(i);
      subplot(121);
      hold on;
      plot(x,y,'mx','LineWidth',2);
      subplot(122);
      hold on;
      plot(x,y,'mx','LineWidth',2);
      ind=(round(x)-1)*s+round(y);
      tgt=eye(4);
      actaxis(i,:)=axis(ind,:);
      tgtaxis(i,:)=axis(ind,:)-paxis*dot(axis(ind,:),paxis);
      tgt(2:end,2:end)=rmat(tgtaxis(i,:),pi/2);
      fidp(i)=trace(results(ind).pm*tgt')/4;
      fprintf('Fidelity of %g,%g is %g\n',x,y,fidp(i));
  end  
    out.figs=[out.figs 221]    
    [out.(name).data.fid ind]= max(fidp);
    % ind is the index into the arrays made above.
    out.(name).bestx = xp(ind);
    out.(name).besty = yp(ind);
    out.(name).tgtaxis = tgtaxis(ind,:);
    out.(name).actaxis = actaxis(ind,:);
    % mind is the index into the main data
    mind=(round(xp(ind))-1)*s+round(yp(ind));
    out.(name).data.chisq = matrix.pdata(mind).chisq;
    out.(name).data.approx_fid = matrix.pdata(mind).approx_fid;
    out.(name).data.eigs = matrix.pdata(mind).eigs;
    subplot(121);
    hold on;
    plot(xp(ind),yp(ind),'mo','LineWidth',3);
    subplot(122);
    hold on;
    plot(xp(ind),yp(ind),'mo','LineWidth',3);
end

%Internal functions below
%%
% Recenter the dot and lock the gradient.
function PrepareMeasurement(opts)
global tunedata;
% Always start by making sure the dot is centered
if isopt(opts,'ana')
    return;
end
pause(1);
awgcntrl('on start wait err raw');
pause(1);
done=0;
while ~done
    autotune('stp tl tmp');
    done=all(tunedata.measp(1:2,1) < .05e-3);
    
    if all(tunedata.measp(1:2,1) < .5e-3)
        atcenter('noconfirm');
    else
        fprintf('Warning: large change.  Will wait for user confirmation\n');
        atcenter;
    end
    pause(1);
end
 if ~sm_setgradient
     error('Unable to lock gradient.  Is fbdata.gradtgt correct?\n');
 end
end

%% Check that the gradient is still locked
function vl=VerifyLock(opts)
if isopt(opts,'ana')
    vl=1;
    return;
end 
  global fbdata;
  g=sm_getgradient([],struct('nloop',100));
  vl = abs(g-fbdata.gradtgt(str2num(opts.datachan(end))));
end

%% Prepare to exit
function out=cleanup(out,opts)
save(out.fname,'out');
if ~isopt(opts,'ana')
  sleep;
end
end

function attr_log(out, logstr, args)
if ~exist('args','var') 
   args = ''; 
end
out.log = [out.log sprintf(logstr,args)];
fprintf(logstr,args);

end

%%
function out=RunProcessTest(out, name, tname, tgt)   
  %now run the max fidelity process over and over to get some stats
  proctime=out.(name).proctime;
  eps = out.(name).data.eps*ones(65,1);
  varpar=[eps, [0 ; repmat(proctime, length(eps)-1,1)]];
  popts=struct('name',tname,'process',struct('type','@exch'),'params',[8 1.15 proctime],'varpar',varpar);
  out=RunProcessTomo(out, opts,popts);
  target = tgt;
  pmat= ana_processTomo_c(out.(tname).fname,struct('tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
  pdesc=ana_process_matrix(pmat, {'opts','eigenbasis', 'tgt',target});
  fids = [pdesc.pdata(2:end).fid];
  out.(tname).fid  = mean(fids);
  out.(tname).data = pdesc;
  out.(tname).vfid = std(fids)/sqrt(length(fids));
  ax=[]; ang=[];
  for h = 1:length(fids)
      ax(:,h) = pdesc.pdata(h+1).eigs.axis;
      ang(h) = pdesc.pdata(h+1).eigs.angle;
  end
  out.(tname).axis = mean(ax,2);
  out.(tname).vaxis = std(ax')/sqrt(length(fids));
  out.(tname).angle = mean(ang);
  out.(tname).vangle = std(ang)/sqrt(length(fids));
  out.(tname).done=1;
end
%% Run process tomo
function out=RunProcessTomo(out, opts, popts)
% popts must have name, process, params, varpar, opts
% popts may have pulsenum
  global awgdata;
    ana = isopt(opts,'ana');
    name=popts.name;
    popts=def(popts,'pulsenum',21856);
    popts=def(popts,'process2','@null');
    popts=def(popts,'process3','@null');
    popts=def(popts,'nrep',25);
    
    if ~ana
        awgrm(opts.lastgroup,'after');
        awgclear('unused');
        % Starting states for process tomography
        preps = struct(... %'ptnull','@null', ...
            'ptad','@adprep', ...
            'ptdephased',struct('type',{'@adprep','@exch'},'time',{[],0.05},'val',{[],[0 0]}), ...
            'ptDU',{{'@dbzpi','@adprep'}}, ...
            'ptpi','@dbzpi', ...
            'pty',{{'@dbzpi','@dbzprep'}}, ...
            'ptdbz','@dbzprep');
        eps100=out.jcal.ffunc_i(out.jcal.fp,100e-3); %epsilon of 100MHz (e-3 is because we use GHz)
        eps60=out.jcal.ffunc_i(out.jcal.fp,60e-3); %epsilon of 60MHz (e-3 is because we use GHz)
        attr_log(out, 'Selected eps=%g for 100MHz J \n',eps100);
        attr_log(out, 'Selected eps=%g for 60MHz J \n',eps60);
        %out.log = [ out.log sprintf('Selected eps=%g for 100MHz J\n',eps100)];
        %out.log = [ out.log sprintf('Selected eps=%g for 60MHz J\n',eps60)];
        
        ptrot=struct('type','@exch','time',0.01,'val',[nan nan eps60]);                
        % Tomographic directions
        tomos={'ST','UD','Y'};
        clear pg;             
        
        % Assemble dictionary entries
        pcal=[];
        gen = sprintf('%03d_a',size(popts.varpar,1));
        pfields=fields(preps);
        for i=1:length(pfields)
            for j=1:length(tomos)
                pg.chan=opts.chan;
                pg.ctrl='pls loop pack';
                pg.dict={opts.side, struct('prep',['@' pfields{i}],'tomo',['@', tomos{j}, 'read'],'process',popts.process,'process2',popts.process2,'process3',popts.process3,'ptrot',ptrot),preps};                
                pg.name=sprintf('attr_processTomo_%s_%d_%d_%s_%s',gen, popts.pulsenum,length(pcal)+1,tomos{j},upper(opts.side(1)));
                pg.pulses=popts.pulsenum; 
                pg.trafofn.func=@skineffect_trafofn;
                pg.trafofn.args=opts.skin;                                               
                % Parameters: p(1) total pulse lenght, p(2) process val(3), p(3) process time(1)
                % n.b. a 2ns pi/2 pulse is about 125 MHz.
                pg.params=popts.params;
                pg.varpar=popts.varpar;
                pcal{end+1}=pg.name;
                
                try
                    plsupdate(pg);
                catch
                    plsdefgrp(pg);
                end
            end
        end
        awgadd(pcal);
        adddbz(upper(opts.side(1)));
        PrepareMeasurement(opts);    
        out.(name).fname=smnext(sprintf(['attr_processTomo_rot_' name '_%s'],upper(opts.side(1))));
        d=smrun(fConfSeq2([length(awgdata.pulsegroups),(opts.lastgroup+1):(length(awgdata.pulsegroups)-1)], ...
            struct('nloop',400,'nrep',popts.nrep,'opts','pol','datachan',opts.datachan)),out.(name).fname);
        out.(name).fname=['sm_' out.(name).fname];
        out.log = [out.log sprintf('  -> Took a %s process Calibration dataset (%s)\n',name, out.tomocal.fname)];
        if any(isnan(d{1}(:)))
            attr_log('    Scan aborted by user\n');
            %out.log=[out.log sprintf('    Scan aborted by user\n')];
            out=cleanup(out,opts);
         	error('abort');
        end
        if ~VerifyLock(opts)
            attr_log('  -> Gradient came unlocked during test!  aborting\n');
            %out.log=[out.log sprintf('  -> Gradient came unlocked during test!  aborting\n')];
            fprintf('Gradient came unlocked during Tomocal test\n');
            out=cleanup(out,opts);
            error('abort');            
        end
    end
end

%% Argument processing
% Apply a default.
function s=def(s,f,v)
if(~isfield(s,f))
    s=setfield(s,f,v);
end
end
function b=isdone(out,name)
  b=isfield(out,name) && isfield(out.(name),'done') && out.(name).done;
end

function b=isopt(config,name)
b=~isempty(strfind(config.opts,name));
end

function b=istest(config,name)
b=~isempty(strfind(config.tests,name));
end

% Add a test requierd by another test
function opts=requires(opts,name,requires)
   if ~iscell(requires)
       requires={requires};
   end
   if ~isempty(strfind(opts.tests,name))
      for i=1:length(requires)
          if isempty(strfind(opts.tests,requires{i}))
              opts.tests=[opts.tests ' ' requires{i}];
          end
      end
   end   
end
