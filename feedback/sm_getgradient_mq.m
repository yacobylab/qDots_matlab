function [grad out]=sm_getgradient(grpname,opts)
% function [grad out]=sm_getgradient(grpname,opts)
%  grpname; name/number of group to use.
%     empty/does not exist: use first group whose name starts with dBz_
%  opts.datachan: default to DAQ2.
%  opts.opts: nopol, nodisp, reget
%  opts.ampthresh: 0.1; If visiblity is less than this, assume fit is bad and
%     make grad_dev huge.
%  opts.chisqthresh: 10; If chi^2 is more than this, assume fit is bad and
%     make grad, grad_dev hughe.

%
%     reget assumes the most recent previous scan to run was sm_getgradient,
%     and skips re-initializing the DAQ card.
% Return the current magnetic field gradient. 
% Sign of the gradient is only trustworthy if the gradient is locked.
% Negative means triplet-side.
% out contains 
%    norm_sig: normalized signal; 1 indicated 100% visibility
%        scan: scan used to take the data
%    grad_dev: estimate of error bar on gradient measurement
%       xd,yd: raw data
%       ff,fp: fit funtion, fit parameters.



%% fill in default options, make the scan
  global fbdata;
  global awgdata;
  if ~exist('opts','var')
      opts=struct();
  end  
  opts=def(opts,'datachan',{'DAQ1','DAQ2'}); 
  opts=def(opts,'chisqthresh',5);
  opts=def(opts,'ampthresh',.2); 
  opts=def(opts,'fitwrapopts','');
  opts=def(opts,'opts','');
  opts=def(opts,'nloop',200);
  if ~iscell(opts.datachan)
      opts.datachan={opts.datachan};
  end
  if ~exist('grpname','var') || isempty(grpname)
     if length(opts.datachan) > 1   
       dbzgrps=find(~cellfun('isempty',regexp({awgdata(1).pulsegroups.name},'^dBz_.*LR')));
     else
       dbzgrps=find(~cellfun('isempty',regexp({awgdata(1).pulsegroups.name},'^dBz_')));
     end
     grpname=dbzgrps(1);  % Default to first DBZ group.
  end
  
  
  if isstruct(grpname)
      scan=grpname;
      grpname=scan.data.pulsegroups(1).name;
  else
      if isopt(opts,'nopol')
        scan=fConfSeq2(grpname,struct('nloop',opts.nloop,'nrep',1,'opts','raw'));
      else
        scan=fConfSeq2(grpname,struct('nloop',opts.nloop,'nrep',1,'opts','pol raw'));
      end
      scan.datachan=opts.datachan; % Protect cell arrays.
  end
  if isopt(opts,'nodisp') 
      scan.disp=[];
  else
      for i=1:length(opts.datachan)
        scan.disp=struct('loop',1,'channel',i,'dim',1);
      end
  end
      
  scan.configch=[];
 
  %% Take the data  
  if isopt(opts,'reget')
    smatrigfn(1,smchaninst(scan.loops(1).getchan{1}),4);
    for i=2:length(scan.loops(1).prefn)
      if ischar(scan.loops(1).prefn(i).fn)
        f=str2func(scan.loops(1).prefn(i).fn);
        f(1,scan.loops(1).prefn(i).args{:});          
      else
        scan.loops(1).prefn(i).fn(1,scan.loops(1).prefn(i).args{:});
      end
    end
    for i=1:length(opts.datachan)
      rd=smget(scan.loops(1).getchan{i});
      smset('PulseLine',1);
      d(i).hist=histc(rd{1},scan.loops(1).procfn(4).fn.args{1})';
      d(i).data=mean(reshape(rd{1},scan.loops(1).procfn(1).fn(3).args{1}),2)';      
      d(i).raw=rd{1};
    end
  else
    rd=smrun(scan);
    nc=length(opts.datachan);
    for i=1:length(opts.datachan)
        d(i).hist=nc+i+1;
        d(i).data=i;
        d(i).raw=nc+i+2;
    end
  end
  if any(isnan(d(1).data))
      grad=nan;
      return;
  end
  
  % Cache the xvals for faster repeated scans.
  if isfield(scan.data,'xv')
      xv=scan.data.xv;
  else
      xv=plsinfo('xval',grpname);
      scan.data.xv = xv;
  end
  
  %% Loop through the channels and do the fitting.
  for chan=1:length(d)
      % We guess the channel number by looking at the last digit of the
      % getchan.  This will not scale well to more than 9 qubits!
      cnum=str2num(opts.datachan{chan}(end))
      % Work out x,y, sigma_y
      data_std=std(reshape(d(i).raw,length(d(i).data),opts.nloop),0,2)/sqrt(opts.nloop);           
      out(chan).scan=scan;
      out(chan).xd=xv;
      out(chan).yd=squeeze(d(chan).data);
      
      fbdata.refval(cnum) = out.yd(1);
      
      % Initial guess for fit
      fp=fioscill(xv,squeeze(d(chan).data),1);
      fp(5:6)=0;
      cosfn2 = '@(y, x)y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*exp(-(x-y(5)).^2 * y(6).^2)';
      
      % Actual fit
      [fp,chisq,cov]=mfitwrap(struct('x',xv,'y',d(chan).data,'vy',(data_std.^2)'),struct('fn',str2func(cosfn2)),fp,'',[1 1 1 1 0 0]);
      
      % Work out the visibility of the oscillation
      vc=scan.loops(1).procfn(end).fn.args{1};
      vc=vc+(vc(2)-vc(1))/2;
      n=sum(d(chan).hist);
      xbar = sum(d(chan).hist .* vc)/n;
      xxbar = sum(d(chan).hist .* vc .* vc)/n;
      sigstd=sqrt(xxbar-xbar*xbar);
      out(chan).norm_sig=sqrt(fp(2)^2+fp(3)^2)/sigstd;
      
      % Use derivative of fit function to guess singlet or triplet side
      dcosfn = @(y,x) -y(2)*sin(y(4)*x) + y(3) * cos(y(4)*x);
      fbtime=fbdata.taus(cnum);
      grad(chan) = 1e3*(fp(4)/(2*pi))*sign(dcosfn(fp,fbtime));
      out(chan).grad_dev = 1e3*sqrt(cov(4,4))/(2*pi);
      out(chan).fp=fp;
      out(chan).ff=str2func(cosfn2);
      
      if ~isnan(opts.ampthresh) && out(chan).norm_sig < opts.ampthresh
          out(chan).grad_dev=1000;
          grad(chan)=0;
      end
      if ~isnan(opts.chisqthresh) && chisq > opts.chisqthresh
          grad(chan)=0;
          out.grad_dev=1000;
      end
  end
end

% Apply a default.
function s=def(s,f,v)
if(~isfield(s,f))
    s=setfield(s,f,v);
end
return;
end

function b=isopt(config,name)
b=~isempty(strfind(config.opts,name));
return;
end