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
%   t1 -- allow t1 decay in process matricezs.
%   nophys -- allow non-physical process matrices.
% opts.lastgroup: TuneRotations will delete groups after this.  Default 23
% opts.side: default "right"
% opts.chan: default [3 4]
% opts.skin: skin effect comp (default 0)
% opts.tests: string containing
%      t1 -- check t1
% tomocal -- calibrates tomography.  Implies t1


    global tunedata;    
    global awgdata;
    
    addpath z:/oliver/matlab/bloch/
    addpath z:/qdots/matlab/tomocal/
        
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
    opts.awgq = awgdata.quiet;
    awgdata.quiet = ~isopt(opts,'chatty');

    opts=def(opts,'tests','tunes');
    opts=def(opts,'side','right');
    opts=def(opts,'chan',[3 4]);       % fixme; get from tunedata
    opts=def(opts,'datachan','DAQ2');  % fixme; get from tunedata
    opts=def(opts,'lastgroup',23);
    opts=def(opts,'skin',0); % Skin effect @ 1 Ghz in dB
    opts=requires(opts,'all',{'tests','tunes','proccal_load'});
    opts=requires(opts,'tests',{'proctest_j_plus','proctest_j_minus','proctest_j_pi'});
    opts=requires(opts,'tunes',{'proccal_j_pi','proccal_hjh','proccal_hjjh',...
        'proccal_j_plus','proccal_j_minus','proccal_jjh','proccal_hjmh','proccal_j4','proccal_hjj',...
        'proccal_hjjhjj','proccal_h'});
    if ~ana % fixme; loop on this until it doesn't change
        opts=requires(opts,'proccal_jjh',{'proccal_j_pi','proccal_hjh'});
        opts=requires(opts,'proccal_h','proccal_hjh');
        opts=requires(opts,'proccal_hjh','proccal_j_plus');
        opts=requires(opts,'proccal_hjjh',{'proccal_j_pi','proccal_hjh'});
        opts=requires(opts,'proccal_hjmh','proccal_j_minus');
        opts=requires(opts,'proctest_j_pi','proccal_j_pi');
        opts=requires(opts,'proccal_j_pi',{'tomocal','proccal_j_plus'});
        opts=requires(opts,'proctest_j_minus','proccal_j_minus');
        opts=requires(opts,'proccal_j_minus',{'tomocal','proccal_j_plus'});
        opts=requires(opts,'proctest_j_plus','proccal_j_plus');
        opts=requires(opts,'proccal_j4','tomocal');
        opts=requires(opts,'proccal_j8','tomocal');
        opts=requires(opts,'proccal_j_plus','tomocal');
        opts=requires(opts,'proccal_b_plus','tomocal');
        opts=requires(opts,'proccal_load','tomocal');
        opts=requires(opts,'tomocal',{'t1','jcal'});
        opts=requires(opts,'jcal','t1');
    end    
    
    if ~exist('out','var') || isempty(out)
        out.log=sprintf('Started a new at_TuneRotations run at %s\n',datestr(now));
        out.fname=smnext('attr_summary');
        out.figs=[];        
    else
        out.log=[out.log sprintf('Continued at_TuneRotations run at %s\n',datestr(now))];
    end
    if ~isfield(out,'opts')
        out.opts=opts;
    end
    out=attr_log(out,'Tests to run: %s\n%g dB skin comp\n',opts.tests,opts.skin);
    
    % If ana is set, don't run any scans or change any gates.
    ana = isopt(opts,'ana');
    if ana
        anastr=' ana';
    else
        anastr='';
    end
    awgrm(opts.lastgroup,'after');
    awgclear('unused');
    
    opts=def(opts,'apm_arg','');
    if isopt(opts,'t1') % options to pass to ana_process_matrix
        opts.apm_arg=[opts.apm_arg ' t1'];
    end
    if isopt(opts,'nophys') % options to pass to ana_process_matrix
        opts.apm_arg=[opts.apm_arg ' nophys'];
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
                out=attr_log(out, '  -> Gradient came unlocked during test!  aborting\n');
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
            out=attr_log(out,'  -> Tunerun %d.%d: measured T1=%g, vis=%g\n',[out.t1.run,out.t1.subrun, out.t1.t1, out.t1.vis']);
            out.t1.done=1;
        end
        out.figs=[out.figs [79 80 81]];
        tmpdict=pdload(tunedata.name);
        out.dict.meas = tmpdict.meas;
        out.dict.meas.time(1) = 1e6*tunedata.runs(out.t1.run).t1(out.t1.subrun).tmeas;
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
            for eps=(.7:.1:2.5)+.1
                pg.xval = eps;
                pg.params=[eps 0];
                pg.trafofn.func=@skineffect_trafofn;
                pg.trafofn.args=opts.skin;
                pg.varpar=(0:1:63)';
                pg.name=sprintf(namepat,length(jgrp)+1);
                try
                    plsupdate(pg);
                catch e
                    fprintf('Error adding group: %s: %s\n',e.identifier,e.message);    
                    plsdefgrp(pg);
                end
                jgrp=[jgrp pg.name];
            end
            fprintf('loading jcal pulses /n');
            awgadd(jgrp);
            adddbz(upper(opts.side(1)));
            
            PrepareMeasurement(opts);
            
            out.jcal.fname=smnext(sprintf('attr_Ramsey_%s',upper(opts.side(1))));
            d=smrun(fConfSeq2([length(awgdata.pulsegroups),(opts.lastgroup+1):(length(awgdata.pulsegroups)-1)], ...
                struct('nloop',400,'nrep',8,'opts','pol','datachan',opts.datachan)),out.jcal.fname);  
            out.jcal.fname=['sm_' out.jcal.fname];
            out=attr_log(out,'  -> Took a Ramsey dataset (%s)\n',out.jcal.fname);
            if any(isnan(d{1}(:)))
                out=attr_log(out,'    Scan aborted by user\n');
                out=cleanup(out,opts);
                return
            end
            if ~VerifyLock(opts)
                out=attr_log(out,'  -> Gradient came unlocked during test!  aborting\n');
                out=cleanup(out,opts);
                return;
            end
            
        end
        [figs fits data]=ana_echo(out.jcal.fname,struct('opts','noppt ramsey freq afitdecay','rng',[6 inf],'fb',300));
        gooddata=find(fits.freq(2,:) > 200 & imag(fits.freq(2,:)) == 0,1,'last'):size(fits.freq,2);   % These are in MHz.
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
            slowfudge=.1;
            fastfudge=0;
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
                        pg.varpar=(1:32)' -1;
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
                    catch e
                        fprintf('Error adding group: %s: %s\n',e.identifier,e.message);    
                        plsdefgrp(pg);
                    end
                    tcalg=[tcalg pg.name];
                end
            end
            fprintf('loading tomocal pulses\n')
            awgadd(tcalg);        
            adddbz(upper(opts.side(1)));            
            PrepareMeasurement(opts);
            out.tomocal.fname=smnext(sprintf('attr_TomoCal_rot_%s',upper(opts.side(1))));
            d=smrun(fConfSeq2([length(awgdata.pulsegroups),(opts.lastgroup+1):(length(awgdata.pulsegroups)-1)], ...
                struct('nloop',400,'nrep',10,'opts','pol','datachan',opts.datachan)),out.tomocal.fname);  
            out.tomocal.fname=['sm_' out.tomocal.fname];
            out=attr_log(out,'  -> Took a Tomography Calibration dataset (%s)\n',out.tomocal.fname);
            if any(isnan(d{1}(:)))
                out=attr_log(out,'    Scan aborted by user\n');
                out=cleanup(out,opts);
                return
            end
            if ~VerifyLock(opts)
                out=attr_log(out,'  -> Gradient came unlocked during test!  aborting\n');
                fprintf('Gradient came unlocked during Tomocal test\n');
                out=cleanup(out,opts);
                return;
            end            
        end
        out.tomocal.cal=ana_tomocal_2(out.tomocal.fname, struct('opts','ref quiet noppt'));
        out.figs=[out.figs out.tomocal.cal.figs];
        out.tomocal.done=1;        
    end        
    
    %% Now tune up %J-rotations with process tomo
    if istest(opts,'proccal_load') && (ana || ~isdone(out,'proccal_load'))
        % params are pulse length, epsilon ,pulse time
        varpar=linspace(0,.1,65)';
        
        popts=struct('name','proccal_load','process',struct('type','@reload'),'params',[nan nan .1],'varpar',varpar,'pulsenum',21857);
        out=RunProcessTomo(out, opts, popts);
        
        dummy  = ana_processTomo_c(out.proccal_load.fname,struct('opts','','tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
        target = [1 0 0 0; 0 0 0 0; 0 0 0 0; 1 0 0 0]; % Load singlet.
        out.proccal_load.pdata = ana_process_matrix(dummy, {'opts',opts.apm_arg, 'tgt',target});
        out.proccal_load.done=1;
    end
    out=def(out,'axes',struct('x',[1 0 0],'y',[0 1 0],'z',[0 0 1]));
    out.axes=def(out.axes,'x',[1 0 0]);
    out.axes=def(out.axes,'y',[0 1 0]);
    out.axes=def(out.axes,'z',[0 0 1]);
                
    %% For our J-rotation, ignore the axis and just get the angle right.
    name='proccal_j_plus';
    if istest(opts,name) && ((ana && isfield(out,name)) || (~ana && ~isdone(out,name)))
        % params are pulse length, epsilon ,pulse time
        proctime=2e-3;
        eps = [0, linspace(out.jcal.ffunc_i(out.jcal.fp,70e-3),...
            out.jcal.ffunc_i(out.jcal.fp,300e-3), 64)]';
        varpar=[eps, [0 ; repmat(proctime, length(eps)-1,1)]];
        
        popts=struct('name',name,'process',struct('type','@exch'),'params',[8 1.15 proctime],'varpar',varpar);
        
        out=RunProcessTomo(out, opts, popts);
        
        % Extract the process matrix
        pmat = ana_processTomo_c(out.(name).fname,struct('opts','quiet','tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
        % Analyze the process into a rotation and dephasing
        pdesc=ana_process_matrix(pmat, {'opts',opts.apm_arg});
        for i=1:length(pdesc.pdata)
            ang(i)=pdesc.pdata(i).rot.angle;
            axis(i,:)=pdesc.pdata(i).rot.axis;
            if axis(i,3) < 0
                axis(i,:) = -axis(i,:);
                ang(i) = -ang(i);
            end
        end
        [~, ind] = min(abs(ang-pi/2));                    % Pick the rotation closest to pi/2
        out.axes.z = axis(ind,:);                     % this *defines* the z axis.
        out.(name).target = rmat4(out.axes.z, pi/2);   % make the target matrix and maximize fidelity 
        pdesc=ana_process_matrix(pmat, {'opts',opts.apm_arg,'tgt',out.(name).target});
        [~, indmaxfid] = max([pdesc.pdata.fid]);
        if indmaxfid ~= ind
            out = attr_log(out,'  +*** Warning: maximum fidelity rotation is unexpected (%d ~= %d)\n',ind, indmaxfid);
        end
        
        out.(name).data = pdesc.pdata(ind);
        out.(name).data.eps = eps(ind);
        out.(name).data.fid = pdesc.pdata(ind).fid;
        out.(name).data.proctime = proctime;
        out.(name).proc = 'R_z \pi/2';
        out.(name).dictentry = 'Rz_pi2';
        out.dict.(out.(name).dictentry) = ...
                   struct('type','@exch','time',out.(name).data.proctime,'val',[nan nan out.(name).data.eps]);
        out.(name).done=1;
        out = attr_log(out,'  + Fidelity of %s = %g\n',name,out.(name).data.fid);
    end
    tname='proctest_j_plus';
    if istest(opts,tname) && (ana || ~isdone(out,tname))
        out=RunProcessTest(out, name, tname);  
    end
    
    name = 'proccal_hjh';
    if istest(opts,name) && ((ana && isfield(out,name)) || (~ana && ~isdone(out,name)))
        % params are pulse length, epsilon ,pulse time
        numjs = 16; %number of Js to try (loop iterations)]
        numjps= 16;
        jrng = .35;
        jprng = .35;
        jtime = out.proccal_j_plus.data.proctime;
        
        % Warning: trig below.  Oliver may suck.
        %theta=(atan2(sqrt(sum(out.axes.z(1:2).^2)), abs(out.axes.z(3))) + .35*pi);
        theta=atan2(sqrt(sum(out.axes.z(1:2).^2)), abs(out.axes.z(3)));
        nu_dbz = out.jcal.fits.omega_dbz * 1e3 / (2*pi);        
        %jpnu = nu_dbz / tan(theta);        
        jpnu = nu_dbz * tan(pi/4-theta);
        if jpnu < 0
            error('Ouch!');
        end
        
        if jpnu-8 < 1e3*out.jcal.fp(3)
            jpnu = 1e3*out.jcal.fp(3)+20; %25MHz more than the min freq            
            out = attr_log(out,['  + Desired frequency is less than offset in jcal.\n',...
                '  +Choosing jpnu=%g MHz \n'],jpnu);
            fprintf('frequency less than jcal offset requested\n');
        end
        jptime = 1e-3*ceil(1e3*1/((2*sqrt(2)*jpnu)));
        out=attr_log(out,'  + dBz=%gh MHz.  Targetting %g MHz %g ns H rotation\n', nu_dbz, jpnu, jptime*1e3);
        
        jval  = -.3+linspace(out.proccal_j_plus.data.eps-.5*jrng, out.proccal_j_plus.data.eps+.5*jrng,numjs);
        jptgt=1e-3*jpnu;
        jpval = linspace(out.jcal.ffunc_i(out.jcal.fp,jptgt)-.4*jprng, out.jcal.ffunc_i(out.jcal.fp,jptgt)+.6*jprng,numjps);
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
        out=attr_ana_jprime(out,name,out.axes.z,opts);
        out.axes.x = out.(name).tgtaxis;
        out.axes.y = cross(out.axes.z, out.axes.x);
        out.(name).data.jeps=interp1(jval,out.(name).bestx);
        out.(name).data.heps=interp1(jpval,out.(name).besty);
        out.(name).data.jtime = jtime(2);
        out.(name).data.htime = jptime(2);        
        out.(name).proctime = jptime(2)*2+jtime(2);
        out.(name).dictentry = 'Rx_pi2';
         out.dict.(out.(name).dictentry) = ...
                  struct('type',{'@exch','@exch','@exch'},...
                          'time',{out.(name).data.htime,out.(name).data.jtime,out.(name).data.htime},...
                          'val',{[nan nan out.(name).data.heps],[nan nan out.(name).data.jeps],[nan nan out.(name).data.heps]});
        out.(name).target=rmat4(out.axes.x, pi/2);
        out.(name).done=1;
        figure(221);
        out = attr_log(out,'  + Fidelity of %s = %g\n',name,out.(name).data.fid);
    end      
    
    name='proccal_h';   target = rmat4((out.axes.z + out.axes.x)/2, pi);
    if istest(opts,name) && ((ana && isfield(out,name)) || (~ana && ~isdone(out,name)))    
        % params are pulse length, epsilon ,pulse time
        proctime=1e-3*floor(out.proccal_hjh.data.htime*1e3/2);
        eps = [linspace(out.proccal_hjh.data.heps-.3,out.proccal_hjh.data.heps+.3, 16)]';
        [eps1 eps2] = meshgrid(eps,eps);
        eps1   = [0; eps1(:)];        
        eps2  = [0; eps2(:)];        
        times  = [0;  repmat(proctime,length(eps1)-1,1)];
        
        varpar = [ eps1, times, eps2, times+1e-3];
        
        popts=struct('name',name,'nrep', 15,...
            'process',struct('type','@exch'),...
            'process2',struct('type','@exch'), ...
            'params',[5 1.15 1e-3 1 1e-3],'varpar',varpar,'pulsenum',21858);        
        
        out=RunProcessTomo(out, opts, popts);
        
        % Extract the process matrix
        pmat = ana_processTomo_c(out.(name).fname,struct('opts','quiet','tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
        % Analyze the process into a rotation and dephasing
        pdesc=ana_process_matrix(pmat, {'opts',opts.apm_arg, 'tgt',target});
        [maxfid indmaxfid]=max([pdesc.pdata.fid]);        
        out.(name).data = pdesc.pdata(indmaxfid);
        out.(name).data.eps1 = eps1(indmaxfid);
        out.(name).data.eps2 = eps2(indmaxfid);
        out.(name).data.time = proctime;
        out.(name).data.fid = maxfid;
        out.(name).data.proctime = proctime*2;
                out.(name).proc = 'H';
        out.(name).dictentry = 'R_H';
          out.dict.(out.(name).dictentry) = ...
                  struct('type',{'@exch','@exch'},'time',{proctime, proctime}, ...
                  'val',{[nan nan out.(name).data.eps1],[nan nan out.(name).data.eps2]});
        out.(name).target=target;
        out.(name).done=1;
        out = attr_log(out,'  + Fidelity of %s = %g\n',name,out.(name).data.fid);
    end
    
    name='proccal_j_pi';   target = rmat4(out.axes.z, pi);
    if istest(opts,name) && ((ana && isfield(out,name)) || (~ana && ~isdone(out,name)))
        % params are pulse length, epsilon ,pulse time
        proctime=3e-3;
        eps = [0, linspace(out.jcal.ffunc_i(out.jcal.fp,70e-3),...
            out.jcal.ffunc_i(out.jcal.fp,300e-3), 64)]';
        varpar=[eps, [0 ; repmat(proctime, length(eps)-1,1)]];
        
        popts=struct('name',name,'process',struct('type','@exch'),'params',[8 1.15 proctime],'varpar',varpar);
        
        out=RunProcessTomo(out, opts, popts);
        
        % Extract the process matrix
        pmat = ana_processTomo_c(out.(name).fname,struct('opts','quiet','tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
        % Analyze the process into a rotation and dephasing
        pdesc=ana_process_matrix(pmat, {'opts',opts.apm_arg, 'tgt',target});
        [maxfid indmaxfid]=max([pdesc.pdata.fid]);        
        out.(name).data = pdesc.pdata(indmaxfid);
        out.(name).data.eps = eps(indmaxfid);
        out.(name).data.fid = maxfid;
        out.(name).data.proctime = proctime;
                out.(name).proc = 'R_z \pi';
        out.(name).dictentry = 'Rz_pi';
          out.dict.(out.(name).dictentry) = ...
                  struct('type','@exch','time',out.(name).data.proctime,'val',[nan nan out.(name).data.eps]);
        out.(name).target=target;
        out.(name).done=1;
        out = attr_log(out,'  + Fidelity of %s = %g\n',name,out.(name).data.fid);
   end
    tname='proctest_j_pi';
    if istest(opts,tname) && (ana || ~isdone(out,tname))
        out=RunProcessTest(out, name, tname);  
    end
    
    name='proccal_j_minus';  target =  rmat4(out.axes.z, 3*pi/2);
    if istest(opts,name) && ((ana && isfield(out,name)) || (~ana && ~isdone(out,name)))
        % params are pulse length, epsilon ,pulse time
        proctime=5e-3;
        eps = [0, linspace(out.jcal.ffunc_i(out.jcal.fp,70e-3),...
            out.jcal.ffunc_i(out.jcal.fp,300e-3), 64)]';
        varpar=[eps, [0 ; repmat(proctime, length(eps)-1,1)]];
        
        popts=struct('name',name,'process',struct('type','@exch'),'params',[8 1.15 proctime],'varpar',varpar);
        out=RunProcessTomo(out, opts, popts);
        
        pmat  = ana_processTomo_c(out.(name).fname,struct('opts','','tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
        pdesc=ana_process_matrix(pmat, {'opts',opts.apm_arg, 'tgt',target});
        [maxfid indmaxfid]=max([pdesc.pdata.fid]);
        out.(name).data = pdesc.pdata(indmaxfid);
        out.(name).data.eps = eps(indmaxfid);
        out.(name).data.proctime=proctime;
        out.(name).data.fid = maxfid;
        out.(name).dictentry = 'Rz_3pi2';
         out.dict.(out.(name).dictentry) = ...
                  struct('type','@exch','time',out.(name).data.proctime,'val',[nan nan out.(name).data.eps]);
        out.(name).target=target;
        out.(name).done=1;
        out = attr_log(out,'  + Fidelity of %s = %g\n',name,out.(name).data.fid);
   end
    tname='proctest_j_minus';
    if istest(opts,tname) && (ana || ~isdone(out,tname))
        out=RunProcessTest(out, name, tname);  
    end            
   
    name='proccal_hjmh';  target=rmat4(out.axes.x, 3*pi/2);
    if istest(opts,name) && ((ana && isfield(out,name)) || (~ana && ~isdone(out,name)))
        % params are pulse length, epsilon ,pulse time
        numjs = 16; %number of Js to try (loop iterations)]
        numjps= 16;
        jrng = .4;
        jprng = .4;
        jtime = out.proccal_j_minus.data.proctime;
        jptime = out.proccal_hjh.data.htime;

        jval  = linspace(out.proccal_j_minus.data.eps-.5*jrng, out.proccal_j_minus.data.eps+.5*jrng,numjs);
        %jpval = linspace(out.jcal.ffunc_i(out.jcal.fp,35e-3)-.5*jprng, out.jcal.ffunc_i(out.jcal.fp,35e-3)+.5*jprng,numjps);
        jpval = -.4+linspace(out.proccal_hjh.data.heps-.5*jprng,out.proccal_hjh.data.heps+.5*jprng,numjps);
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
        pmat  = ana_processTomo_c(out.(name).fname,struct('opts','','tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
        pdesc=ana_process_matrix(pmat, {'opts',opts.apm_arg, 'tgt',target});

        [maxfid indmaxfid]=max([pdesc.pdata.fid]);
        out.(name).data = pdesc.pdata(indmaxfid);
        out.(name).data.jeps = epsj(indmaxfid);
        out.(name).data.heps = epsjp(indmaxfid);
        out.(name).data.jtime = jtime(2);
        out.(name).data.htime = jptime(2);
        out.(name).data.proctime=2*jptime(2)+jtime(2);
        out.(name).data.fid = maxfid;
        out.(name).dictentry = 'Rx_3pi2';
         out.dict.(out.(name).dictentry) = ...
           struct('type',{'@exch','@exch','@exch'},...
                   'time',{out.(name).data.htime,out.(name).data.jtime,out.(name).data.htime},...
                   'val',{[nan nan out.(name).data.heps],[nan nan out.(name).data.jeps],[nan nan out.(name).data.heps]});
        out.(name).target=target;
        out.(name).done=1;
         out = attr_log(out,'  + Fidelity of %s = %g\n',name,out.(name).data.fid);
    end
    
    name='proccal_hjjh';  target=rmat4(out.axes.x, pi);
    if istest(opts,name) && ((ana && isfield(out,name)) || (~ana && ~isdone(out,name)))
        % params are pulse length, epsilon ,pulse time
        numjs = 16; %number of Js to try (loop iterations)]
        numjps= 16;
        jrng = .3;
        jprng = .3;
        jtime = out.proccal_j_pi.data.proctime;
        % For stability it's beneficial to make J' too close to Z axis.
        jptime = out.proccal_hjh.data.htime;

        jval  = -.2+linspace(out.proccal_j_pi.data.eps-.5*jrng, out.proccal_j_pi.data.eps+.5*jrng,numjs);
        %jpval = linspace(out.jcal.ffunc_i(out.jcal.fp,35e-3)-.5*jprng, out.jcal.ffunc_i(out.jcal.fp,35e-3)+.5*jprng,numjps);
        jpval = .1+linspace(out.proccal_hjh.data.heps-.5*jprng,out.proccal_hjh.data.heps+.5*jprng,numjps);
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
        pmat = ana_processTomo_c(out.(name).fname,struct('opts','','tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
        pdesc=ana_process_matrix(pmat, {'opts',opts.apm_arg, 'tgt',target});

        [maxfid indmaxfid]=max([pdesc.pdata.fid]);
        out.(name).data = pdesc.pdata(indmaxfid);
        out.(name).data.jeps = epsj(indmaxfid);
        out.(name).data.heps = epsjp(indmaxfid);
        out.(name).data.jtime = jtime(2);
        out.(name).data.htime = jptime(2);
        out.(name).data.proctime=2*jptime(2)+jtime(2);
        out.(name).data.fid = maxfid;
        out.(name).dictentry = 'Rx_pi';
        out.dict.(out.(name).dictentry) = ...
            struct('type',{'@exch','@exch','@exch'},...
                   'time',{out.(name).data.htime,out.(name).data.jtime,out.(name).data.htime},...
                   'val',{[nan nan out.(name).data.heps],[nan nan out.(name).data.jeps],[nan nan out.(name).data.heps]});
        out.(name).target=target;
        out.(name).done=1;
        out = attr_log(out,'  + Fidelity of %s = %g\n',name,out.(name).data.fid);
    end    
    
    name='proccal_jjh';  target=rmat4(out.axes.y,pi/2); % j-pi followed by H- should give Y pi/2
    if istest(opts,name) && ((ana && isfield(out,name)) || (~ana && ~isdone(out,name)))
        % params are pulse length, epsilon ,pulse time
        numjs = 16; %number of Js to try (loop iterations)]
        numjps= 16;
        jrng = .3;
        jprng = .3;
        jtime = out.proccal_j_pi.data.proctime;
        % For stability it's beneficial to make J' too close to Z axis.
        jptime = out.proccal_hjh.data.htime;

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
        
        pmat  = ana_processTomo_c(out.(name).fname,struct('opts','','tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
        pdesc=ana_process_matrix(pmat, {'opts',opts.apm_arg, 'tgt',target});
        
        
        [maxfid indmaxfid]=max([pdesc.pdata.fid]);
        out.(name).data = pdesc.pdata(indmaxfid);
        out.(name).data.jeps = epsj(indmaxfid);
        out.(name).data.heps = epsjp(indmaxfid);
        out.(name).data.jtime = jtime(2);
        out.(name).data.htime = jptime(2);
        out.(name).data.proctime=2*jptime(2)+jtime(2);
        out.(name).data.fid = maxfid;
        out.(name).dictentry = 'Ry_pi2';
        out.dict.(out.(name).dictentry) = ...
            struct('type',{'@exch','@exch'},...
                   'time',{out.(name).data.jtime,out.(name).data.htime},...
                   'val',{[nan nan out.(name).data.jeps],[nan nan out.(name).data.heps]});
        out.(name).target=target;
        out.(name).done=1;
        out = attr_log(out,'  + Fidelity of %s = %g\n',name,out.(name).data.fid);
    end    
       
    name='proccal_hjj';  target=rmat4(out.axes.y, -pi/2); % j-pi followed by H- should give Y pi/2
    if istest(opts,name) && ((ana && isfield(out,name)) || (~ana && ~isdone(out,name)))
        % params are pulse length, epsilon ,pulse time
        numjs = 16; %number of Js to try (loop iterations)]
        numjps= 16;
        jrng = .3;
        jprng = .3;
        jtime = out.proccal_j_pi.data.proctime;
        % For stability it's beneficial to make J' too close to Z axis.
        jptime = out.proccal_hjh.data.htime;

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
       
        pmat  = ana_processTomo_c(out.(name).fname,struct('opts','','tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
        pdesc=ana_process_matrix(pmat, {'opts',opts.apm_arg, 'tgt',target});
        
        [maxfid indmaxfid]=max([pdesc.pdata.fid]);
        out.(name).data = pdesc.pdata(indmaxfid);
        out.(name).data.jeps = epsj(indmaxfid);
        out.(name).data.heps = epsjp(indmaxfid);
        out.(name).data.jtime = jtime(2);
        out.(name).data.htime = jptime(2);
        out.(name).data.proctime=2*jptime(2)+jtime(2);
        out.(name).data.fid = maxfid;
        out.(name).dictentry = 'Ry_3pi2';
        out.dict.(out.(name).dictentry) = ...
            struct('type',{'@exch','@exch'},...
                   'time',{out.(name).data.htime,out.(name).data.jtime},...
                   'val',{[nan nan out.(name).data.heps],[nan nan out.(name).data.jeps]});
        out.(name).target=target;
        out.(name).done=1;
        out = attr_log(out,'  + Fidelity of %s = %g\n',name,out.(name).data.fid);
    end    
    
    name='proccal_hjjhjj';  target=rmat4(out.axes.y, pi); % j-pi followed by H- should give Y pi/2
    if istest(opts,name) && ((ana && isfield(out,name)) || (~ana && ~isdone(out,name)))
        % params are pulse length, epsilon ,pulse time
        numjs = 16; %number of Js to try (loop iterations)]
        numjps= 16;
        jrng = .3;
        jprng = .3;
        jtime = out.proccal_j_pi.data.proctime;
        % For stability it's beneficial to make J' too close to Z axis.
        jptime = out.proccal_hjh.data.htime;

        jval  = linspace(out.proccal_j_pi.data.eps-.5*jrng, out.proccal_j_pi.data.eps+.5*jrng,numjs);
        jpval = linspace(out.jcal.ffunc_i(out.jcal.fp,35e-3)-.5*jprng, out.jcal.ffunc_i(out.jcal.fp,35e-3)+.5*jprng,numjps);
        [j jp]=meshgrid(jval,jpval);
        epsj   = [0; j(:)];        
        epsjp  = [0; jp(:)];        
        jtime  = [0;  repmat(jtime,length(epsj)-1,1)];
        jptime = [0; repmat(jptime,length(epsj)-1,1)];
        varpar = [ epsjp, jptime, epsj, jtime, epsjp, jptime, epsj, jtime  ];
            
        popts=struct('name',name,'nrep', 15,...
            'process',struct('type','@exch'),...
            'process2',struct('type','@exch'), ...
            'process3',struct('type','@exch'), ...
            'process4',struct('type','@exch'), ...
            'params',[5 1.15 1e-3 1 1e-3 1 1e-3 1 1e-3 ],'varpar',varpar,'pulsenum',21861);
        out=RunProcessTomo(out, opts, popts);
        
        pmat  = ana_processTomo_c(out.(name).fname,struct('opts','','tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
        pdesc=ana_process_matrix(pmat, {'opts',opts.apm_arg, 'tgt',target});
        
        [maxfid indmaxfid]=max([pdesc.pdata.fid]);
        out.(name).data = pdesc.pdata(indmaxfid);
        out.(name).data.jeps = epsj(indmaxfid);
        out.(name).data.heps = epsjp(indmaxfid);
        out.(name).data.jtime = jtime(2);
        out.(name).data.htime = jptime(2);
        out.(name).data.proctime=2*jptime(2)+jtime(2);
        out.(name).data.fid = maxfid;
        out.(name).dictentry = 'Ry_pi';
        out.dict.(out.(name).dictentry) = ...
            struct('type',{'@exch','@exch','@exch','@exch'},...
                   'time',{out.(name).data.htime,out.(name).data.jtime,out.(name).data.htime,out.(name).data.jtime},...
                   'val',{[nan nan out.(name).data.heps],[nan nan out.(name).data.jeps],[nan nan out.(name).data.heps],[nan nan out.(name).data.jeps]});
        out.(name).target=target;
        out.(name).done=1;
        out = attr_log(out,'  + Fidelity of %s = %g\n',name,out.(name).data.fid);
    end        
    
    name='proccal_j4';  target =  rmat4(out.axes.z, pi/4);
    if istest(opts,name) && (ana || ~isdone(out,name))
        % params are pulse length, epsilon ,pulse time
        proctime=1e-3;
        eps = [0, linspace(out.jcal.ffunc_i(out.jcal.fp,70e-3),...
            out.jcal.ffunc_i(out.jcal.fp,300e-3), 64)]';
        varpar=[eps, [0 ; repmat(proctime, length(eps)-1,1)]];
        
        popts=struct('name',name,'process',struct('type','@exch'),'params',[8 1.15 proctime],'varpar',varpar);
        out=RunProcessTomo(out, opts, popts);
        
        pmat  = ana_processTomo_c(out.(name).fname,struct('opts','','tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
        pdesc=ana_process_matrix(pmat, {'opts',opts.apm_arg, 'tgt',target});
        [maxfid indmaxfid]=max([pdesc.pdata.fid]);
        out.(name).data = pdesc.pdata(indmaxfid);
        out.(name).data.eps = eps(indmaxfid);
        out.(name).data.proctime=proctime;
        out.(name).data.fid = maxfid;
        out.(name).dictentry = 'Rz_pi4';
         out.dict.(out.(name).dictentry) = ...
                  struct('type','@exch','time',out.(name).data.proctime,'val',[nan nan out.(name).data.eps]);
        out.(name).target=target;
        out.(name).done=1;
        out = attr_log(out,'  + Fidelity of %s = %g\n',name,out.(name).data.fid);
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
function out = attr_ana_jprime(out,name,paxis,opts)
% function out = attr_ana_jprime(out,name,paxis)
%  paxis is axis to be perpendicular to.
  file=out.(name).fname;
  tomofile=[out.tomocal.fname '_mat.mat'];
  
  process=ana_processTomo_c(file,struct('opts','','tomocal',tomofile));
  matrix=ana_process_matrix(process,{'opts',opts.apm_arg});
  
  results=matrix.pdata(2:end);
  for i=1:length(results)
      angle(i)=results(i).rot.angle;
      axis(i,:)=results(i).rot.axis;
      dp(i)=dot(results(i).rot.axis,paxis); 
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
  cmap_center; colormap(cmap_str('yr.gb')); colorbar;
  subplot(122);
  imagesc(dp); hold on; colorbar; title('Axis Error');
  cmap_center; colormap(cmap_str('yr.gb')); colorbar;
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
      actaxis(i,:)=axis(ind,:);
      tgtaxis(i,:)=axis(ind,:)-paxis*dot(axis(ind,:),paxis);
      tgt=rmat4(tgtaxis(i,:),pi/2);
      fidp(i)=trace(results(ind).pm*tgt')/4;
      fprintf('Fidelity of %g,%g is %g\n',x,y,fidp(i));
  end  
    out.figs=[out.figs 221];   
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
    out.(name).data.rot = matrix.pdata(mind).rot;
    out.(name).data.pm = matrix.pdata(mind).pm;
    subplot(121);
    hold on;
    plot(xp(ind),yp(ind),'mo','LineWidth',3);
    subplot(122);
    hold on;
    plot(xp(ind),yp(ind),'mo','LineWidth',3);
end
% Internal functions below.

%%
% Recenter the dot and lock the gradient.
function PrepareMeasurement(opts)
if isopt(opts,'ana')
    return;
end
attr_PrepareMeasurement();
end
%% Check that the gradient is still locked
function vl=VerifyLock(opts)
if isopt(opts,'ana')
    vl=1;
    return;
end 
  global fbdata;
  g=sm_getgradient([],struct('nloop',100));
  vl = abs(g-fbdata.gradtgt(str2double(opts.datachan(end))));
end

%% Prepare to exit
function out=cleanup(out,opts)
save(out.fname,'out');
if ~isopt(opts,'ana')
  sleep;
end
global awgdata;
awgdata.quiet = opts.awgq;
end

function out=attr_log(out, logstr, varargin)
out.log = [out.log sprintf(logstr,varargin{:})];
fprintf(logstr,varargin{:});
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
  pmat= ana_processTomo_c(out.(tname).fname,struct('opts','quiet','tomocal',[out.tomocal.fname '_mat.mat'])); %extract the process matrix
  pdesc=ana_process_matrix(pmat, {'opts',opts.apm_arg, 'tgt',target});
  fids = [pdesc.pdata(2:end).fid];
  out.(tname).fid  = mean(fids);
  out.(tname).data = pdesc;
  out.(tname).vfid = std(fids)/sqrt(length(fids));
  ax=[]; ang=[];
  for h = 1:length(fids)
      ax(:,h) = pdesc.pdata(h+1).rot.axis;
      ang(h) = pdesc.pdata(h+1).rot.angle;
  end
  out.(tname).axis = mean(ax,2);
  out.(tname).vaxis = std(ax')/sqrt(length(fids));
  out.(tname).angle = mean(ang);
  out.(tname).vangle = std(ang)/sqrt(length(fids));
  out.(tname).done=1;
  save(out.fname,'out');
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
    popts=def(popts,'process4','@null');

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
        
        ptrot=struct('type','@exch','time',0.01,'val',[nan nan eps60]);                
        if ~isfield(out,'dict') || ~isfield(out.dict,'ptrot')
            out.dict.ptrot = ptrot;
        end
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
                pg.dict={opts.side, struct('prep',['@' pfields{i}],'tomo',['@', tomos{j}, 'read'],'process',popts.process,'process2',popts.process2,'process3',popts.process3,'process4',popts.process4,'ptrot',ptrot),preps};                
                pg.name=sprintf('attr_processTomo_g2_%s_%d_%d_%s_%s',gen, popts.pulsenum,length(pcal)+1,tomos{j},upper(opts.side(1)));
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
                catch e
                    fprintf('Error adding group: %s: %s\n',e.identifier,e.message);    
                    plsdefgrp(pg);
                end
            end
        end
        fprintf('adding  groups for %s\n',popts.name);
        awgadd(pcal);
        adddbz(upper(opts.side(1)));
        PrepareMeasurement(opts);    
        out.(name).fname=smnext(sprintf(['attr_processTomo_rot_' name '_%s'],upper(opts.side(1))));
        d=smrun(fConfSeq2([length(awgdata.pulsegroups),(opts.lastgroup+1):(length(awgdata.pulsegroups)-1)], ...
            struct('nloop',400,'nrep',popts.nrep,'opts','pol','datachan',opts.datachan)),out.(name).fname);
        out.(name).fname=['sm_' out.(name).fname];
        out.log = [out.log sprintf('  -> Took a %s process Calibration dataset (%s)\n',name, out.tomocal.fname)];
        if any(isnan(d{1}(:)))
            out=attr_log('    Scan aborted by user\n');
            %out.log=[out.log sprintf('    Scan aborted by user\n')];
            out=cleanup(out,opts);
         	error('abort');
        end
        if ~VerifyLock(opts)
            out=attr_log('  -> Gradient came unlocked during test!  aborting\n');
            fprintf('Gradient came unlocked during Tomocal test\n');
            out=cleanup(out,opts);
            error('abort');            
        end
    end
    save(out.fname,'out');
end

%% Argument processing
% Apply a default.
function s=def(s,f,v)
if(~isfield(s,f))
    s.(f) = v;
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
