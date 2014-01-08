function out=at_TestProcess(opts,attr)
% function out=at_TestProcess(opts,attr)
% Test a process.
%   opts.name -- name of process to test
%   opts.target (optional) -- name of target process.
%   opts.nloop -- easy.
%   opts.nrep -- number of reps
%   opts.chop -- divide reps into this many subgroups for variances.
%   opts.dict -- dictionary entries that define process.
%   attr: result from running at_tunerotations
%   opts.opts ana

global awgdata;

opts=def(opts,'nrep',20);
opts=def(opts,'nloop',500);
opts=def(opts,'chop',5);
opts=def(opts,'opts','');
opts=def(opts,'length',4);
opts=def(opts,'skin',0); % skin effect compensation
opts=def(opts,'name','Jane_Doe');
name=sprintf('attp_%s',opts.name);

if ~isopt(opts,'ana')
    
    if ~isreal(attr.jcal.fp)
       error('Unphysical Jcal. Consider using a different attr');
    end
    eps60=attr.jcal.ffunc_i(attr.jcal.fp,60e-3); %epsilon of 60MHz (e-3 is because we use GHz)
    eps100=attr.jcal.ffunc_i(attr.jcal.fp,100e-3); %epsilon of 60MHz (e-3 is because we use GHz)
    if ~isreal(eps60) || ~isreal(eps100)
       error('imaginary epsilons requested. Bad idea'); 
    end
    % Get the prep dict ready
    preps = struct(... %'ptnull','@null', ...
        'ptad','@adprep', ...
        'ptdephased',struct('type',{'@adprep','@exch'},'time',{[],0.05},'val',{[],[0 0]}), ...
        'ptDU',{{'@dbzpi','@adprep'}}, ...
        'ptdbz1',struct('type','@sep','time',0.004), ...
        'ptdbz2',struct('type','@sep','time',0.008), ...
        'ptdbz3',struct('type','@sep','time',0.012), ...
        'ptdbz4',struct('type','@sep','time',0.016), ...
        'ptdbz5',struct('type','@sep','time',0.020), ...
        'ptdbz6',struct('type','@sep','time',0.024), ...
        'ptdbz7',struct('type','@sep','time',0.028), ...        
        'ptrot1',struct('type',{'@adprep','@exch'},'time',{[],0.002},'val',{[],[nan nan eps100]}), ...
        'ptrot2',struct('type',{'@adprep','@exch'},'time',{[],0.004},'val',{[],[nan nan eps100]}), ...
        'ptrot3',struct('type',{'@adprep','@exch'},'time',{[],0.006},'val',{[],[nan nan eps100]}), ...
        'ptrot4',struct('type',{'@adprep','@exch'},'time',{[],0.008},'val',{[],[nan nan eps100]}));
        
    f=fieldnames(preps);
    for i=1:length(f)
        preps.(sprintf('prep%02d',i))=['@' f{i}];
    end
    
    out.prepdict = preps;
    out.attr_fname=attr.fname;
    out.sumname=smnext(name);
    out.fname=smnext(sprintf('proctest_%s_%s',opts.name,upper(attr.opts.side(1))));    
    
    if ~iscell(opts.dict)
       for i=1:length(opts.dict(:))
          d{i}=opts.dict(i);
       end
       opts.dict=d;
    end
    opts.dict=[{struct('process','@null')}, opts.dict]; % Add a null process to start of dict for a reference.
    out.dict=opts.dict;
    out.opts = opts;
    out.attr_dict=attr.dict;
    tomos={'ST','UD','Y'};
    allpulses=21884:21904;
    pcal=[];
    for i=1:length(out.dict)
        for j=1:length(tomos)
            pg.chan=attr.opts.chan;
            pg.ctrl='pls loop pack';
            pg.dict={attr.opts.side, out.prepdict,struct('tomo',['@', tomos{j}, 'read'],'ptrot',attr.dict.ptrot), attr.dict, out.dict{i}};
            pg.name=sprintf('attr_processTest_B_%d_%d_%d_%s_%s',opts.length*10,i,length(f),tomos{j},upper(attr.opts.side(1)));
            pg.pulses=allpulses(1:length(f));
            pg.params=[opts.length nan nan];
            pg.trafofn.func=@skineffect_trafofn;
            pg.trafofn.args=opts.skin;
            pg.xval=i-1;
            pcal{end+1}=pg.name;
            try
                plsupdate(pg);
            catch e
                fprintf('Error adding group: %s: %s\n',e.identifier,e.message);
                plsdefgrp(pg);
            end
        end
    end
    awgrm(attr.opts.lastgroup,'after');
    awgclear('unused');    
    awgadd(pcal);
    attr_PrepareMeasurement();
    tmp=out.fname;
    out.fname=['sm_' out.fname];
    save(out.sumname,'out');
    d=smrun(fConfSeq2((attr.opts.lastgroup+1):(length(awgdata.pulsegroups)), ...
            struct('nloop',opts.nloop,'nrep',opts.nrep,'opts','pol','datachan',attr.opts.datachan)),tmp);
    sleep;
end
end

% Argument processing
% Apply a default.
function s=def(s,f,v)
if(~isfield(s,f))
    s.(f) = v;
end
end

function b=isopt(config,name)
b=~isempty(strfind(config.opts,name));
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
